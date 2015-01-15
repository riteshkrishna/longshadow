% 
% sam_file            = 'example-data/Dog-150-autmate/cfam-150-unsplicedHuman.sam'
% bed_sam_file        = 'example-data/Dog-150-autmate/cfam-150-unsplicedHuman-new.bed.fa'
% sam_read_file       = 'example-data/Dog-150-autmate/cfam-150-unsplicedHuman-samread.fa'
% cigar_threshold     = 160
% out_small_sam_txt   = 'example-data/Dog-150-autmate/small-Only-Cfam-150.txt'
% out_cigar_mat       = 'example-data/Dog-150-autmate/dog_150_cigar.mat'
% out_jc69_mat        = 'example-data/Dog-150-autmate/dog_150_jc69.mat'
% out_titv_mat        = 'example-data/Dog-150-autmate/dog_150_titv.mat'
% out_multivar_file   = 'example-data/Dog-150-autmate/dog_150_multivariateData.mat'
% workspace_dir       = 'example-data/Dog-150-autmate'

% Execute pipeline
% mapping_pipeline_run(sam_file,bed_sam_file,sam_read_file,...
%                 cigar_threshold, out_small_sam_txt, out_cigar_mat, ...
%                 out_jc69_mat,out_titv_mat,out_multivar_file,workspace_dir)
%     

sam_file            = 'example-data/Dog-150-autmate/cfam-150-unsplicedHuman.sam'
bed_sam_file        = 'example-data/Dog-150-autmate/cfam-150-unsplicedHuman-new.bed.fa'
out_small_sam_txt   = 'example-data/Dog-150-autmate/small-Only-Cfam-150.txt'
out_data_file       = 'example-data/Dog-150-autmate/results.txt'

%perform_grep_awk_on_sam(sam_file, out_small_sam_txt)
fprintf(1, 'Done Grep')
        
CIGAR_JC69_titv_multivariateData(out_small_sam_txt,bed_sam_file,out_data_file)