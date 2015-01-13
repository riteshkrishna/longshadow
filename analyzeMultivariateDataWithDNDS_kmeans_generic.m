
%% Prepare data for further analysis
    multivariate_mat = 'example-data/Dog-150-autmate/dog_150_multivariateData.mat';
    map_class = 3; 
   [filtered_data,filtered_reads] = analyzeMultivariateData_generic(out_multivar_file,map_class);
        
%% Compute KaKs values for chosen reads
    sam_read_fasta = 'example-data/Dog-150-autmate/cfam-150-unsplicedHuman-samread.fa';
    bed_read_fasta = 'example-data/Dog-150-autmate/cfam-150-unsplicedHuman-new.bed.fa';
       
    [final_result_nan,final_result_dnds] = call_compute_dnds_for_reads_generic(filtered_data,filtered_reads,sam_read_fasta, bed_read_fasta)

    dnds_ids = final_result_dnds.read_ids;
    dnds_values = final_result_dnds.read_data;

    
%% Perform k-means
    X = [dnds_values(:,1)  dnds_values(:,3)]
    idx_right = perform_kmeans(X);
    final_idx1 = dnds_values(idx_right==1,:);
    final_idx2 = dnds_values(idx_right==2,:);
    figure; plotmatrix(final_idx1);
    figure; plotmatrix(final_idx2);
    

    