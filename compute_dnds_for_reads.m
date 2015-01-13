function [reads_kaks_map,sam_reads_seqs_map,bed_reads_seqs_map] = ...
            compute_dnds_for_reads(readlist,sam_read_fasta, bed_read_fasta)

    sam_reads_seqs_map = containers.Map();
    bed_reads_seqs_map = containers.Map();
    reads_kaks_map = containers.Map();

    fid = fopen(sam_read_fasta);
    tline = fgetl(fid);
    while ischar(tline) % Read file one line at a time
        C = strread(tline,'%s','delimiter','\t');
        read_id     = C{1};
        sam_seq     = C{2};
    
        sam_reads_seqs_map(read_id) = sam_seq;
        tline = fgetl(fid);
    end
    fclose(fid);

    % Read read-id and mapped-bed-sequences
    fid = fopen(bed_read_fasta);
    tline = fgetl(fid);
    while ischar(tline) % Read file one line at a time
        C = strread(tline,'%s','delimiter','\t');
        read_id     = C{1};
        sam_seq     = C{2};
    
        bed_reads_seqs_map(read_id) = sam_seq;
        tline = fgetl(fid);
    end
    fclose(fid);
    
    % Compute DNDS for all reads from the readlist
    for i=1:size(readlist)
        this_read = readlist{i};
        sam_seq = sam_reads_seqs_map(this_read);
        bed_seq = bed_reads_seqs_map(this_read);
        KaKs = compute_dnds(sam_seq,bed_seq);
        reads_kaks_map(this_read) = KaKs;
    end

end