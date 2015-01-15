%% This code reads a SAM file and produces a TSV text file with three  columns
% The three columns are - QNAME, RNAME and CIGAR
% The mappings with CIGAR = '*' are ignored while producing the text file.
%
%
% The example of equivalent unix commmand is -
% !grep -v '@' example-data/Dog-150/cfam-150-unsplicedHuman.sam | awk '{ if ( match($1,"Cfam")) {print $0}}' >> example-data/Dog-150//Only-Cfam-150.sam
% ! grep -v '@' example-data/Dog-150//Only-Cfam-150.sam | awk '{if($6 !="*") {print $0}}' | cut -f1,3,6  >> example-data/Dog-150/small-Only-Cfam-150.txt
%
% Example :
% sam_file = 'example-data/Dog-150/Only-Cfam-150.sam'
% out_small_sam_txt = '/Users/ritesh/tmp/readsam.txt'
% perform_grep_awk_on_sam(sam_file, out_small_sam_txt )

%%
function perform_grep_awk_on_sam(sam_file, out_small_sam_txt )

    f_sam = fopen(sam_file, 'r');
    f_out = fopen(out_small_sam_txt,'w');
    
    header_regx = '@([a-zA-Z1-9_]+)';
     
    count = 1;
    tline = fgetl(f_sam);
    while ischar(tline)
        ret = regexp(tline, header_regx);
        if length(ret) > 0
            %fprintf(1,'%s',tline)
            %continue;
        else
            C = strread(tline,'%s','delimiter','\t');
    
            QNAME   = C{1};
            FLAG    = C{2};
            RNAME   = C{3};
            POS     = C{4};
            MAPQ    = C{5};
            CIGAR   = C{6};
            SEQ     = C{10};
            
            TF = strfind(CIGAR,'*');
            if isempty(TF) 
                fprintf(f_out,'%s\t%s\t%s\t%s\n',QNAME, RNAME, CIGAR, SEQ);
            end
            
        end
        
        tline = fgetl(f_sam);
    end

    fclose(f_sam)
    fclose(f_out)

end

