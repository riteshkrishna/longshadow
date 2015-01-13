%% Compute KaKs/DnDs for selected reads
% Fit this prior to run
% sam_read_fasta = '/Users/ritesh/Ritesh_CGR_Work/Ortholog-data/Multi-Species-Approach/FromEnsemble/April/Analysis-Bowtie-JC69/Sample-Sams/DN-DS/example-data/Dog-150/cfam-150-unsplicedHuman-samread.fa';
% bed_read_fasta = '/Users/ritesh/Ritesh_CGR_Work/Ortholog-data/Multi-Species-Approach/FromEnsemble/April/Analysis-Bowtie-JC69/Sample-Sams/DN-DS/example-data/Dog-150/cfam-150-unsplicedHuman-new.bed.fa';

function [final_result_nan,final_result] = call_compute_dnds_for_reads_generic(read_data,read_ids,sam_read_fasta, bed_read_fasta)


    %% Computing begins
    kaks_col = zeros(size(read_ids,1),1);

    kaks_map = compute_dnds_for_reads(read_ids,sam_read_fasta, bed_read_fasta);
    
    for i=1:size(read_ids,1)
        this_read = read_ids{i};
        kaks_col(i) = kaks_map(this_read);
    end

    %% Add the data as the fourth column
    read_data = [read_data kaks_col];

    %% Find out which indices are NaNs
    nan_idx = find(any(isnan(kaks_col),2) == 1);
    inf_idx = find(any(isinf(kaks_col),2) ==1) % Inf indices will have CIGAR =1
    nan_inf_idx = unique([nan_idx;inf_idx]);


    read_data_nan = read_data(nan_inf_idx,:);
    read_ids_nan = read_ids(nan_inf_idx);
    read_ids(nan_inf_idx) = [];
    read_data(nan_inf_idx,:) = [];


    field_1 = 'read_data_nan';
    field_2 = 'read_ids_nan';
    field_3 = 'read_ids';
    field_4 = 'read_data';

    % different order of storing...?
    final_result_nan = struct(field_1,read_data_nan,field_2,read_ids_nan);
    final_result = struct(field_3,read_ids,field_4,read_data);
end
