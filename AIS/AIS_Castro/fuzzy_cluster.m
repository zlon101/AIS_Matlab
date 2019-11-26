function [F] = fuzzy_cluster(M,mCe);

%
% Ph.D. Thesis
% Copyright by Leandro Nunes de Castro
% June, 2000
% Immune Network (iNet) - Description in aiNet.doc
% Fuzzy Clustering from Hard Clustering (centers of mass)
% 
% function [F] = fuzzy_cluster(M,mCe);
% F   -> matrix of fuzzy membership to each cluster (c x N)
% M   -> matrix of memory cells
% mCe -> matrix of centroids
%

[Nc,l] = size(mCe); F = [];
for i=1:Nc,
   F(:,i) = dist(M,mCe(i,:)');
end;
F = norma(1./F');
F = logsig(10.*F);


% Function normalizes matrix over [0,1]
function [Dn] = norma(D);
% Dn  -> normalized vector over [0,1]
[np,ni] = size(D);
if ni == 1,
   Dn = (D - min(D))./(max(D)-min(D));
else,
   vmaxD = max(D); vminD = min(D);
   for i=1:ni,
      Dn(:,i) = (D(:,i) - vminD(i))./(vmaxD(i)-vminD(i));
   end;
end;
% End Function NORMA

