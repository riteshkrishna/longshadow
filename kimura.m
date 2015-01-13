%% Copied from seqpdist.m and modified..

function [d,R]= kimura(s1,s2,indelMethod)
% indelMethod = s, p, c (for score, pairwise-delete, and complelte-delete resp.)
%
% Kimura distance
%
% Let P and Q be the transition and transversional difference proportions:
% d = - log(1 -2*P -Q)/2 - log(1 -2*Q)/4;
%
% Reference: Kimura, M. (1980)  Journal of Molecular Evolution

[vsc,val] = nwalign(s1,s2,'alphabet','NT');
s1 = val(1,:);
s2 = val(3,:);

if indelMethod == 's'  % score gaps
        h = (s1~='-') | (s2~='-'); % only avoid sites with double gaps
 else                   % pairwise deleted gaps
        h = (s1~='-') & (s2~='-'); % avoid any site with gaps
end

L_temp = [nt2int(s1(h),'ACGTOnly', true)' , nt2int(s2(h),'ACGTOnly', true)'];
L_temp(find(ismember(L_temp,0)),:) = [];
L = double(L_temp)

X  = accumarray(L,1,[4 4]);
numPairs = sum(X(:));
numTransitions = sum(X([3 8 9 14])); % 3:G->A, 8:T->C, 9:A->G, 14:C->T
numTransversions = numPairs - sum(diag(X)) - numTransitions;
P = numTransitions / numPairs;
Q = numTransversions / numPairs;
d = - log(1 -2*P -Q)/2 - log(1 -2*Q)/4;

R = numTransitions/numTransversions;
end