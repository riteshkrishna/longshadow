%% Parse a CIGAR string to find the total count of letter of interest. 
% The letters can be M, I, D or S. 
%
% Example : for a CIGAR 30M10I10D20M and letter = M, we should get the
% answer as 50 (as 30 + 20 = 50), the rex_exp = [IDS]
%
% Eg -
% ex_cigar = 4M4I6M7I23M1I8M1D1M1D3M2I3M3I5M
% total_match = parse_cigar_for_total_M_count(ex_cigar,'M','[IDS]')


function [total_match] = parse_cigar_for_total_M_count(cigar_str,letter,reg_exp)
    
    X = strread(cigar_str,'%s','delimiter',letter);
    total_match = 0;
    
    for i=1:size(X)
        m_count = 0;
        
        char_arry = char(X(i));
        idx = regexp(char_arry,reg_exp);
        
        if isempty(idx)
            m_count = str2num(char_arry(1:end));
        else 
            m_count = str2num(char_arry(idx(end)+1:end));
        end
        
        if isempty(m_count)
            total_match = total_match + 0;
        else
            total_match = total_match + m_count;
        end
        
    end
    
end