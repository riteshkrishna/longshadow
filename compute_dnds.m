%%http://www.mathworks.co.uk/products/bioinfo/code-examples.html?file=/products/demos/shipping/bioinfo/dndsdemo.html
function [KaKs] = compute_dnds(sam_seq,bed_seq)
    [score,alignment]= nwalign(sam_seq,bed_seq);
    %seq1 = seqinsertgaps(sam_seq,alignment(1,:))
    %seq2 = seqinsertgaps(bed_seq,alignment(3,:))
    %[dn,ds] = dnds(seq1,seq2)
    [dn,ds] = dnds(sam_seq,bed_seq);
    KaKs=dn/ds
end