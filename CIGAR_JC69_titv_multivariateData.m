% This is an effort to combine code from CIGAR_calculation.m,
% JC69_calculation.m, titv_kimura_calculation.m and
% createMultivariateDataForMat_faster_generic.m

function CIGAR_JC69_titv_multivariateData(sam_file,bed_sam_file,out_data_file)
    
    
    % Data file to contain numbers produced by different methods
    data_storage_file = out_data_file;
    f_out = fopen(data_storage_file,'w');
    
    
    N_blocks = 1e8; % load these many first records from beginning in memory
    bed_block = choose_bed_block(bed_sam_file, N_blocks)
    
    % Read SAM derived file, read-id-mappedRefernces, reads-cigar, read-seq
    fid = fopen(sam_file);
    tline = fgetl(fid);
    line = 0;
    while ischar(tline) % Read file one line at a time
        line = line + 1
        
        C = strread(tline,'%s','delimiter','\t');
        sam_read_id     = C{1};
        ref_source  = C{2};
        cigar       = C{3};
        sam_seq     = C{4};
    
        bed_seq = find_in_bed_block(bed_block,sam_read_id);
        % if not found in above block, then search throughout the file
        if size(bed_seq,1) == 0
            bed_seq = find_in_bed_file(bed_sam_file,sam_read_id);
        end
        
            
        if size(bed_seq,1) ~= 0
            [jc69_distance, titv_ratio] = compute_jc69_titv(sam_seq,bed_seq);
            
            if ~isnan(jc69_distance) & ~isinf(jc69_distance) & ~isnan(titv_ratio) & ~isinf(titv_ratio)
                cigar_score = compute_cigar(cigar);
                %fprintf(f_out,'%s %f %f %f %s %s\n',sam_read_id,cigar_score,titv_ratio,jc69_distance,sam_seq,bed_seq);
                fprintf(f_out,'%s %f %f %f\n',sam_read_id,cigar_score,titv_ratio,jc69_distance);
            end
            
        end 
        tline = fgetl(fid);
    end
    
    fclose(fid);
    fclose(f_out);

end

function bed_block = choose_bed_block(bed_sam_file, N_blocks)
    fileID = fopen(bed_sam_file);
    formatSpec = '%s %s';
    bed_block = textscan(fileID,formatSpec,N_blocks,'Delimiter','\t');
    fclose(fileID);
    
end

function [bed_seq] = find_in_bed_block(bed_block,sam_read_id)
    
    index = find(~cellfun(@isempty,strfind(bed_block{1}, sam_read_id)));
    % if multiple occurances or no occurance, then read doesn't exist
    if (isempty(index) | (numel(index) > 1))
        bed_seq = '';
    else
        bed_seq = bed_block{2}{index};
    end
end

function [bed_seq] = find_in_bed_file(bed_sam_file,sam_read_id)

    bed_seq = '';
    
    % Memory based approach
    N_blocks = 1e8;
    fileID = fopen(bed_sam_file);
    formatSpec = '%s %s';
    while ~feof(fileID)
        C = textscan(fileID,formatSpec,N_blocks,'Delimiter','\t');
        
        index = find(~cellfun(@isempty,strfind(C{1}, sam_read_id)));
        % if multiple occurances or no occurance, then read doesn't exist
        if (isempty(index) | (numel(index) > 1))
            continue;
        else
            bed_seq = C{2}{index};
        end
    end
    
    fclose(fileID)
    
end

function [score] = compute_cigar(cigar_string)
    letter_to_count = 'M';
    reg_ex_cigar = '[IDS]';
    total_M = parse_cigar_for_total_M_count(cigar_string,letter_to_count,reg_ex_cigar);
    
    letter_to_count = 'I';
    reg_ex_cigar = '[MDS]';
    total_I = parse_cigar_for_total_M_count(cigar_string,letter_to_count,reg_ex_cigar);
    
    letter_to_count = 'D';
    reg_ex_cigar = '[IMS]';
    total_D = parse_cigar_for_total_M_count(cigar_string,letter_to_count,reg_ex_cigar);
    
    % Convert them into percentage (0 for no hit, 1 for all hits)
    cigar_length = total_M + total_I + total_D;
    score = total_M / cigar_length;
end

function [jc69_distance, titv_ratio] = compute_jc69_titv(sam_seq,bed_seq)
    seqs = {sam_seq,bed_seq};
    
    method = 'Jukes-Cantor';
    jc69_distance = seqpdist(seqs,'Method','Jukes-Cantor','Alphabet', 'NT', 'PairwiseAlignment', true);
    
    indelMethod = 'p';
    [d,titv_ratio]= kimura(sam_seq,bed_seq,indelMethod);
    
end
