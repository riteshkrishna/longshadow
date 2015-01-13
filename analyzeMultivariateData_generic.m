%% Run createMultivariateDataForMat_faster_generic.m to prepare the base dataset for
%
% Date - May 7, 2014

%% Mail call
function  [filtered_data,filtered_reads] = analyzeMultivariateData_generic(multivariate_mat,whichClass)

    global CIGAR_COL TITV_COL JC69_COL    
    
    CIGAR_COL = 1;
    TITV_COL  = 2;
    JC69_COL  = 3;
    
    multivariate_data = load(multivariate_mat);
    read_data = multivariate_data.combined_data;
    read_ids = multivariate_data.combined_read_ids;
    
    [good_filtered_read_data,good_filtered_read_ids] = select_good_reads(read_data,read_ids);
    [bad_filtered_read_data,bad_filtered_read_ids]= select_bad_reads(read_data,read_ids);
    [confuse_filtered_read_data,confuse_filtered_read_ids]= select_confusing_reads(read_data,read_ids);
    [poor_filtered_read_data,poor_filtered_read_ids] = select_poor_reads(read_data,read_ids);
    
    switch whichClass
        case 1
            filtered_data = good_filtered_read_data;
            filtered_reads = good_filtered_read_ids;
        case 2
            filtered_data = bad_filtered_read_data;
            filtered_reads = bad_filtered_read_ids;   
        case 3
            filtered_data = confuse_filtered_read_data;
            filtered_reads = confuse_filtered_read_ids;   
        case 4
            filtered_data = poor_filtered_read_data;
            filtered_reads = poor_filtered_read_ids;      
        otherwise
            disp('Wrong parameter for whichClass')
    end
end

%% GOOD READS - High CIGAR, High TiTv
function [good_filtered_read_data,good_filtered_read_ids] = select_good_reads(read_data,read_ids)
    
    cigar_min_threshold = 0.6;
    cigar_max_threshold = 2;
    titv_min_threshold  = 1;
    titv_max_threshold = 20;
    
    [good_filtered_read_data,good_filtered_read_ids] ...
        = filter_goof_reads(read_data,read_ids,cigar_min_threshold,cigar_max_threshold,...
        titv_min_threshold,titv_max_threshold);  
    
    str_count = 'GOOD READS : %f < CIGAR < %f, %d < TiTv < %d \n Total selected = %d';
    selected = size(good_filtered_read_data,1);
    sprintf(str_count,cigar_min_threshold,cigar_max_threshold,titv_min_threshold, titv_max_threshold,  selected)
end

%% BAD READS - Low CIGAR, Low TiTv
function [bad_filtered_read_data,bad_filtered_read_ids]= select_bad_reads(read_data,read_ids)
    
    cigar_min_threshold = 0;
    cigar_max_threshold = 0.6;
    titv_min_threshold  = 0;
    titv_max_threshold = 1;
    
    [bad_filtered_read_data,bad_filtered_read_ids] ...
        = filter_goof_reads(read_data,read_ids,cigar_min_threshold,cigar_max_threshold,...
        titv_min_threshold,titv_max_threshold);
    
    str_count = 'BAD READS : %f < CIGAR < %f, %d < TiTv < %d \n Total Selected = %d' ;
    selected = size(bad_filtered_read_data,1);
    sprintf(str_count,cigar_min_threshold,cigar_max_threshold,titv_min_threshold, titv_max_threshold,  selected)
    
end

%% NEITHER GOOD NOT BAD - High CIGAR, Low TiTv
function [confuse_filtered_read_data,confuse_filtered_read_ids]= select_confusing_reads(read_data,read_ids)
    
    cigar_min_threshold = 0.6;
    cigar_max_threshold = 2;
    titv_min_threshold  = 0;
    titv_max_threshold = 1;
    
    [confuse_filtered_read_data,confuse_filtered_read_ids] ...
        = filter_goof_reads(read_data,read_ids,cigar_min_threshold,cigar_max_threshold,...
        titv_min_threshold,titv_max_threshold);
    
    str_count = 'NEITHER GOOD NOT BAD : %f < CIGAR < %f, %d < TiTv < %d \n Total Selected = %d ';
    selected = size(confuse_filtered_read_data,1);
    sprintf(str_count,cigar_min_threshold,cigar_max_threshold,titv_min_threshold, titv_max_threshold,  selected)
    
end

%% POOR READS- Low CIGAR, High TiTv
function [poor_filtered_read_data,poor_filtered_read_ids] = select_poor_reads(read_data,read_ids)
    
    cigar_min_threshold = 0;
    cigar_max_threshold = 0.6;
    titv_min_threshold  = 1;
    titv_max_threshold = 10;
    
    [poor_filtered_read_data,poor_filtered_read_ids] ...
        = filter_goof_reads(read_data,read_ids,cigar_min_threshold,cigar_max_threshold,...
        titv_min_threshold,titv_max_threshold);
    
    str_count = 'POOR READS : %f < CIGAR < %f, %d < TiTv < %d \n Total Selected = %d ';
    selected = size(poor_filtered_read_data,1);
    sprintf(str_count,cigar_min_threshold,cigar_max_threshold,titv_min_threshold, titv_max_threshold,  selected)
    
end


%% Filter reads based on CIGAR and TiTv threshold criteria
function [titv_cigar_filtered_read_data,titv_cigar_filtered_read_ids] ...
        = filter_goof_reads(read_data,read_ids,cigar_min_threshold,cigar_max_threshold,...
        titv_min_threshold,titv_max_threshold)
    
    global CIGAR_COL TITV_COL JC69_COL
    
    [cigar_filtered_read_data,cigar_filtered_read_ids] = filter_basedOnThreshold(read_data,read_ids,cigar_min_threshold,cigar_max_threshold,CIGAR_COL);
    [titv_cigar_filtered_read_data,titv_cigar_filtered_read_ids] = filter_basedOnThreshold(cigar_filtered_read_data,cigar_filtered_read_ids,titv_min_threshold,titv_max_threshold,TITV_COL);
end

%% Filtering based on threshold >= min && threshold < max for a given column
function [filtered_read_data,filtered_read_ids] ...
                = filter_basedOnThreshold(read_data,read_ids,min_threshold,max_threshold,column)
                
    rows_with_matchedCriteria = find((read_data(:,column) >= min_threshold & read_data(:,column) < max_threshold) == 1);
    right_col_1 = read_data(rows_with_matchedCriteria,1);
    right_col_2 = read_data(rows_with_matchedCriteria,2);
    right_col_3 = read_data(rows_with_matchedCriteria,3);
    filtered_read_data = [right_col_1 right_col_2 right_col_3];
    filtered_read_ids = read_ids(rows_with_matchedCriteria);
end


