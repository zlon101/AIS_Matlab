function T = cluster(Z, cutoff, depth)
%CLUSTER Construct clusters from LINKAGE output.
%   T = CLUSTER(Z, CUTOFF) constructs clusters from cluster tree Z.  Z
%   is a matrix of size M-1 by 3, generated by LINKAGE. CUTOFF is a
%   threshold for cutting the hierarchical tree generated by LINKAGE
%   into clusters. When 0 < CUTOFF < 1, clusters are formed when
%   inconsistent values are greater than CUTOFF (see INCONSISTENT).
%   When CUTOFF is an integer and CUTOFF >= 1, then CUTOFF is
%   considered as the maximum number of clusters to keep in the
%   hierachical tree generated by LINKAGE.
%
%   T = CLUSTER(Z, CUTOFF, DEPTH) construct clusters from cluster
%   tree Z, with inconsistent value CUTOFF as the threshhold and
%   DEPTH as the depth in the inconsitent value computation.
%
%   The output T is a vector of size M that contains cluster number
%   for each observation in the original data.
%
%   See also PDIST, LINKAGE, COPHENET, INCONSISTENT, CLUSTERDATA.
    
%   ZP You, 3-10-98
%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 1.4 $

if nargin < 3
   depth = 2;
end

if nargin < 2
   error('Not enough input arguments.');
end

if (cutoff <= 0) | ((cutoff > 1) & fix(cutoff) ~= cutoff)
   error('The second argument is incorrect.');
end

m = size(Z,1)+1;
T = zeros(m,1);

if cutoff < 1 % inconsistency cutoff
   Y = inconsistent(Z,depth);
   Z(:,4) = Y(:,4); % we only need Y(:,4), and Z(:,1:2). Z(:,3) will be used for storing.
   Z = checkcut(Z, cutoff, m-1);
   T = labeltree(Z, T, m-1, 1);
else % maximum number of clusters cutoff
   if m <= cutoff
      T = (1:m)';
   else
      clsnum = 1;
      for k = (m-cutoff+1):(m-1)
         i = Z(k,1); % left tree
         if i <= m % original node, no leafs
            T(i) = clsnum;
            clsnum = clsnum + 1;
         elseif i < (2*m-cutoff+1) % created before cutoff, search down the tree
            T = clusternum(Z, T, clsnum, i-m);
            clsnum = clsnum + 1;
         end
         i = Z(k,2); % right tree
         if i <= m  % original node, no leafs
            T(i) = clsnum;
            clsnum = clsnum + 1;
         elseif i < (2*m-cutoff+1) % created before cutoff, search down the tree
            T = clusternum(Z, T, clsnum, i-m);
            clsnum = clsnum + 1;
         end
      end
   end
end

function T = clusternum(X, T, c, k)
% assign leaves under cluster c to c.

m = size(X,1)+1;
i = X(k,1); % left tree
if i <= m % original node, no leafs
   T(i) = c;
else % created before cutoff, search down the tree
   T = clusternum(X, T, c, i-m);
end
i = X(k,2); % right tree
if i <= m % original node, no leafs
   T(i) = c;
else % created before cutoff, search down the tree
   T = clusternum(X, T, c, i-m);
end

% helper function, call itself recursively.
function [T, c] = labeltree(X, T, k, c)
% assign leaf node to a cluster number

n = size(X,1)+1;
i = X(k,1);
j = X(k,2);

% left tree first.
if i > n % node i is not a leave, it has subtrees 
   [T, c] = labeltree(X, T, i-n, c); % trace back left subtree 
else
   T(i) = c;
end

% if cut off at k, increase cluster number by 1.
if X(k,3) == 0
   c = c+1;
end

if j > n % node i is not a leave, it has subtrees 
   [T, c] = labeltree(X, T, j-n, c); % trace back left subtree 
else
   T(j) = c;
end


function [X, d] = checkcut(X, s, k)

n = size(X,1)+1;
i = X(k,1); % left subtree index 
j = X(k,2); % right subtree index
if X(k,4) > s
   a = 0;
else
   a = 1;
end

if i > n
   [X, b] = checkcut(X,s,i-n);
else
   b = 1;
end

if j > n
   [X, c] = checkcut(X,s,j-n);
else
   c = 1;
end

d = a*b*c;
X(k,3) = d;
