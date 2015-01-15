% The program creates a multivariate data-structure by combining results
% from individual mat files obtained by computing CIGAR, JC69 and TiTv
% values. This datastructure can be used for further processing.

% Date - May 7, 2014
%
% jc69_MatFile = 'example-data/Dog-150/dog_150_jc69.mat';
% cigar_MatFile = 'example-data/Dog-150/dog_150_cigar.mat';
% titv_MatFile = 'example-data/Dog-150/dog_150_titv.mat';
% out_mat_file = 'example-data/Dog-150/dog_150_multivariateData.mat';

function  [nooutlier_combined_data,out_mat_file] ...
        = createMultivariateDataForMat_faster_generic(jc69_MatFile,cigar_MatFile,titv_MatFile,out_mat_file)


    cigar_data = load(cigar_MatFile);      
    cigar_score_map = cigar_data.cigar_score_reads_map;
    
    titv_data = load(titv_MatFile); 
    titv_score_map = titv_data.titv_reads_map;
    
    jc69_data = load(jc69_MatFile); 
    jc69_score_map = jc69_data.jc69_reads_map;
    
    
    % In some cases, it may not be possible to comoute the JC69 or TiTv
    % values for some of the reads, so in those cases, the computed values
    % will be missing. However, there will always be a CIGAR score for each
    % mapping, so we can use cigar_container as reference for creating a
    % combined dataset and drop those reads where the JC69 or TiTv values
    % are found missing.
    reads_to_consider = size(cigar_score_map,1)
    reads_pool = cigar_score_map.keys;
    
    combined_read_ids = cell(reads_to_consider,1);
    combined_data = zeros(reads_to_consider,3);
    
    missing_reads = {};
    for i=1:size(reads_pool,2)
        read_id = reads_pool{i};
        
        found = 0;
        if (isKey(cigar_score_map,read_id) ~= 0) 
            if (isKey(titv_score_map,read_id) ~= 0) 
                if (isKey(jc69_score_map,read_id) ~= 0) 
                    cigar_score_this_read = cigar_score_map(read_id);
                    titv_score_this_read = titv_score_map(read_id);
                    jc69_score_this_read = jc69_score_map(read_id);
                    
                    found = 1;
                    combined_read_ids{i} = read_id;
                    combined_data(i,:) = [cigar_score_this_read titv_score_this_read jc69_score_this_read];
                end
            end
        end
        
        if(found == 0)
            missing_reads{end+1} = read_id;
        end
    end
    
    %figure; hist( cell2mat(cigar_score_map.values(missing_reads)))
    
    % There will be certain read-ids that will not be present in all the
    % three maps, so, there will be some empty rows for certain indices
    % above. We will need to remove those indices otherwise they will
    % create problems during later function calls.
    bad_idx = find(cellfun(@ischar,combined_read_ids) == 0);
    combined_read_ids(bad_idx) = [];
    combined_data(bad_idx,:) = [];
    
    % remove Inf, NaN from right-mapped dataset
    nan_idx = find(any(isnan(combined_data),2) == 1)
    inf_idx = find(any(isinf(combined_data),2) ==1)
    nan_inf_idx = unique([nan_idx;inf_idx]);
    combined_data(nan_inf_idx,:) = [];
    combined_read_ids(nan_inf_idx) = []; 
    
    % Outlier removals
    [nooutlier_combined_data,I] = fun_removeOutliers(combined_data);
    right_idx = find(all(I,2));
    nooutlier_combined_read_ids = combined_read_ids(right_idx);
    
    %save(out_mat_file,'combined_read_ids','combined_data','nooutlier_combined_data','nooutlier_combined_read_ids');
    savefast(out_mat_file,'combined_read_ids','combined_data','nooutlier_combined_data','nooutlier_combined_read_ids');
end
