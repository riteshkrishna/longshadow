%% Same like DNDS_calculation.m, but counts M, I and D from CIGAR and solely
% depends on the information in the SAM file, and doesn't use the BED file.

function CIGAR_calculation(sam_file,cigar_threshold, out_mat_file)

%sam_file = 'example-data/small-cfam-custom_4_not-a-option.txt';
%cigar_threshold = 80;
%number_of_reads_to_process = 50000;

% The following needs to be given in a different way, not possible to use
% '$' in the following command
%! grep -v '@' $sam_file | awk '{if($6 !="*") {print $0}}' | cut -f1,3,6  >> $extract_sam_file

% Run the cut command
%! grep -v '@' example-data/test.sam | awk '{if($6 !="*") {print $0}}' | cut -f1,3,6  >> example-data/small-test.txt

sam_reads_cigar_map  = containers.Map();
sam_reads_refs_map  = containers.Map();

% Read SAM derived file, read-id-mappedRefernces, reads-cigar
fid = fopen(sam_file);
tline = fgetl(fid);
while ischar(tline) % Read file one line at a time
    C = strread(tline,'%s','delimiter','\t');
    
    read_id     = C{1};
    ref_source  = C{2};
    cigar       = C{3};
    
    sam_reads_cigar_map(read_id) = cigar;
    sam_reads_refs_map(read_id) = ref_source;
    
    tline = fgetl(fid);
end
fclose(fid);


% Select the reads where CIGAR indicates of larger difference between the
% reads and the mappings
read_list_to_investigate = containers.Map();
all_reads = keys(sam_reads_cigar_map);
letter_to_count = 'M';
reg_ex_cigar = '[IDS]';

for i=1:size(all_reads,2)
%for i=1:number_of_reads_to_process                        % Test for a small set
    this_read = all_reads{i};
    this_cigar = sam_reads_cigar_map(this_read);
    total_M = parse_cigar_for_total_M_count(this_cigar,letter_to_count,reg_ex_cigar);
    
    if total_M < cigar_threshold
        read_list_to_investigate(this_read) = this_cigar;
    end
end

% Iterate through read_list_to_investigate and check for the JC69
% calculations
suitable_reads = keys(read_list_to_investigate);
cigar_score_reads_map  = containers.Map();
for i=1:size(suitable_reads,2)
    read_id = suitable_reads{i};
    
    this_cigar = sam_reads_cigar_map(read_id);
    
    letter_to_count = 'M';
    reg_ex_cigar = '[IDS]';
    total_M = parse_cigar_for_total_M_count(this_cigar,letter_to_count,reg_ex_cigar);
    
    letter_to_count = 'I';
    reg_ex_cigar = '[MDS]';
    total_I = parse_cigar_for_total_M_count(this_cigar,letter_to_count,reg_ex_cigar);
    
    letter_to_count = 'D';
    reg_ex_cigar = '[IMS]';
    total_D = parse_cigar_for_total_M_count(this_cigar,letter_to_count,reg_ex_cigar);
    
    % Convert them into percentage (0 for no hit, 1 for all hits)
    cigar_length = total_M + total_I + total_D;
    score = total_M / cigar_length;
    
    cigar_score_reads_map(read_id) = score;
    
end

% Save the last bit
%save('cigar_score-calc.mat', 'cigar_score_reads_map', 'read_list_to_investigate','sam_reads_cigar_map','sam_reads_refs_map')
%save(out_mat_file, 'cigar_score_reads_map', 'read_list_to_investigate','sam_reads_cigar_map','sam_reads_refs_map')
%save(out_mat_file, 'cigar_score_reads_map')
savefast(out_mat_file, 'cigar_score_reads_map');

%[readnames_for_rightwrong_matrix,score_rightwrong_matrix,...
%    right_values,wrong_values,errorCheck_reads_refs_map] = postprocessing_right_wrong_identification(cigar_score_reads_map,sam_reads_refs_map)
%
%plot_rightWrongvalues(score_rightwrong_matrix,right_values,wrong_values)

end