%% Transition/Transversion calculation for read-sequences and mapped-sequences
% Also computes Kimura distance.

% parameters -
% sam_file = SAM file
% bed_sam_file = BED file corresponding to the SAM file 
% sam_read_file = SAM read specific file to go with the BED file
% cigar_threshold = CIGAR threshold to use
% out_mat_file = matlab file to store results
%
% @date - modified on Feb 18, 2014


function titv_kimura_calculation(sam_file,bed_sam_file, sam_read_file, cigar_threshold, out_mat_file)

%sam_file = 'example-data/small-cfam-custom_4_not-a-option.sam';
%bed_sam_file = 'example-data/cfam-custom_4_not-a-option-new-out.fa';
%sam_read_file = 'example-data/cfam-custom_4_not-a-option-new.fa';
%cigar_threshold = 80;
%number_of_reads_to_process = 50000;

% The following needs to be given in a different way, not possible to use
% '$' in the following command
%! grep -v '@' $sam_file | awk '{if($6 !="*") {print $0}}' | cut -f1,3,6  >> $extract_sam_file

% Run the cut command
%! grep -v '@' example-data/test.sam | awk '{if($6 !="*") {print $0}}' | cut -f1,3,6  >> example-data/small-test.txt

sam_reads_seqs_map = containers.Map();
bed_reads_seqs_map = containers.Map();
sam_reads_cigar_map  = containers.Map();
sam_reads_refs_map  = containers.Map();
kimura_reads_map = containers.Map();    % contains Kimura distance
titv_reads_map = containers.Map();      % contains transition/transversion

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

% Read read-id and mapped-bed-sequences
fid = fopen(bed_sam_file);
tline = fgetl(fid);
while ischar(tline) % Read file one line at a time
    C = strread(tline,'%s','delimiter','\t');
    
    read_id     = C{1};
    seq  = C{2};
    
    bed_reads_seqs_map(read_id) = seq;
    tline = fgetl(fid);
end
fclose(fid);

% Read read-id and read-sequences 
fid = fopen(sam_read_file);
tline = fgetl(fid);
while ischar(tline) % Read file one line at a time
    C = strread(tline,'%s','delimiter','\t');
    
    read_id = C{1};
    seq = C{2};
    
    sam_reads_seqs_map(read_id) = seq;
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

% Iterate through read_list_to_investigate and check for the kimura
% calculations
suitable_reads = keys(read_list_to_investigate);
for i=1:size(suitable_reads,2)
    read_id = suitable_reads{i};
    
    if isKey(bed_reads_seqs_map,read_id) == 0
        kimura_reads_map(read_id) = NaN;
        titv_reads_map(read_id) = NaN;
        continue;
    end
    
    read_seq = sam_reads_seqs_map(read_id)         % The read-sequence
    mapped_ref_seq = bed_reads_seqs_map(read_id)   % The ref sequence
    
    %seqs = {read_seq,mapped_ref_seq};
    %SeqsMultiAligned = multialign(seqs,'terminalGapAdjust',true)
    
    indelMethod = 'p';
    [d,R]= kimura(read_seq,mapped_ref_seq,indelMethod)
    
    kimura_reads_map(read_id) = d;
    titv_reads_map(read_id) = R;
end

% Save the last bit
%save(out_mat_file, 'kimura_reads_map', 'titv_reads_map', 'read_list_to_investigate','sam_reads_cigar_map','sam_reads_seqs_map','sam_reads_refs_map')
%save(out_mat_file, 'titv_reads_map')
savefast(out_mat_file, 'titv_reads_map');

% Analysis for Kimura
%[readnames_for_rightwrong_matrix_kimora,score_rightwrong_matrix_kimora,...
%    right_values_kimora,wrong_values_kimora,errorCheck_reads_refs_map_kimora] = postprocessing_right_wrong_identification(kimura_reads_map,sam_reads_refs_map)

%plot_rightWrongvalues(score_rightwrong_matrix_kimora,right_values_kimora,wrong_values_kimora)

% Analysis for Ti-Tv
%[readnames_for_rightwrong_matrix,score_rightwrong_matrix,...
%    right_values,wrong_values,errorCheck_reads_refs_map] = postprocessing_right_wrong_identification(titv_reads_map,sam_reads_refs_map)
%
%plot_rightWrongvalues(score_rightwrong_matrix,right_values,wrong_values)

end

