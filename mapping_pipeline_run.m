%% Runs all the stages in the mapping one by one.
% The parameters can be set in mapping_pipeline_parameters.m file

function mapping_pipeline_run(sam_file,bed_sam_file,sam_read_file,...
                cigar_threshold, out_small_sam_txt, out_cigar_mat, ...
                out_jc69_mat,out_titv_mat,out_multivar_file,workspace_dir)
        
        % Read SAM and extract relevant columns 
        
        perform_grep_awk_on_sam(sam_file, out_small_sam_txt)
        fprintf(1, 'Done Grep')
        
        % Compute CIGAR
        CIGAR_calculation(out_small_sam_txt,cigar_threshold, out_cigar_mat)
        fprintf(1, 'Done CIGAR')
        
        % Compute JC69
        JC69_calculation(out_small_sam_txt,bed_sam_file, sam_read_file, cigar_threshold, out_jc69_mat);
        fprintf(1, 'Done JC69')
        
        % Compute TiTv
        titv_kimura_calculation(out_small_sam_txt,bed_sam_file, sam_read_file, cigar_threshold, out_titv_mat);
        fprintf(1, 'Done TiTv')   
        
        % Create multivariate dataset
        %createMultivariateDataForMat_faster(out_jc69_mat,out_cigar_mat,out_titv_mat,out_multivar_file)
        createMultivariateDataForMat_faster_generic(out_jc69_mat,out_cigar_mat,out_titv_mat,out_multivar_file)
        fprintf(1, 'Done Multivariate Grouping')
        
        % Analyse  multivariate dataset
        map_class = 3
        %[right_values, wrong_values, right_ids, wrong_ids] = analyzeMultivariateData(out_multivar_file, map_class)
        [filtered_data,filtered_reads] = analyzeMultivariateData_generic(out_multivar_file,map_class);
        filtered_data_mat = strcat(workspace_dir,'/Filtered_Data.mat');
        %save(filtered_data_mat,'filtered_data','filtered_reads');
        savefast(filtered_data_mat,'filtered_data','filtered_reads');
        fprintf(1, 'Finished filtering of reads for DNDS')
        
        % Compute KaKs values for chosen reads
        sam_read_fasta = sam_read_file;
        bed_read_fasta = bed_sam_file;
        
        [final_result_nan,final_result_dnds] = call_compute_dnds_for_reads_generic(filtered_data,filtered_reads,sam_read_fasta, bed_read_fasta);
        nan_dnds_ids = final_result_nan.read_ids_nan;
        nan_dnds_values = final_result_nan.read_data_nan;
        dnds_ids = final_result_dnds.read_ids;
        dnds_values = final_result_dnds.read_data;
        
        filtered_dnds_mat = strcat(workspace_dir,'/Filtered_DNDS.mat');
        %save(filtered_dnds_mat,'nan_dnds_ids','nan_dnds_values','dnds_ids','dnds_values');
        savefast(filtered_dnds_mat,'nan_dnds_ids','nan_dnds_values','dnds_ids','dnds_values');
        fprintf(1, 'Finished DNDS')
        
        % Perform k-means
        X = [dnds_values(:,1)  dnds_values(:,3)]
        
        opts = statset('Display','final');
        [idx,ctrs] = kmeans(X,2,'Distance','city','Replicates',5,'Options',opts);
    
        figure;
        plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)

        hold on
        plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
        plot(ctrs(:,1),ctrs(:,2),'kx','MarkerSize',12,'LineWidth',2)
        plot(ctrs(:,1),ctrs(:,2),'ko','MarkerSize',12,'LineWidth',2)
        legend('Cluster 1','Cluster 2','Centroids','Location','NW')
        hold off
        kmeans_figure = strcat(workspace_dir,'/kmeans_fig.eps');
        print('-dtiff','-r300',kmeans_figure)
        
        final_idx1 = dnds_values(idx==1,:);
        final_idx2 = dnds_values(idx==2,:);
        
        figure; plotmatrix(final_idx1);
        cluster1_figure = strcat(workspace_dir,'/cluster1_fig.eps');
        print('-dtiff','-r300',cluster1_figure)
        
        figure; plotmatrix(final_idx2);
        cluster2_figure = strcat(workspace_dir,'/cluster2_fig.eps');
        print('-dtiff','-r300',cluster2_figure)
        
        kmeans_mat = strcat(workspace_dir,'/kmeans.mat');
        %save(kmeans_mat,'X','final_idx1','final_idx2');
        savefast(kmeans_mat,'X','final_idx1','final_idx2');
        fprintf(1, 'Finished K-mean')
        
end

