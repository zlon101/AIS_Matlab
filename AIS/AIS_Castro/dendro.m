function [Z,H,Cn] = dendro(M);

%
% Ph.D. Thesis
% Leandro Nunes de Castro
% February, 2000
% Immune Network (iNet) - Description in iNet.doc
% Secondary Functions: PDIST, LINKAGE, INCONSISTENT, DENDROGRAM
% The secondary functions were taken from the STATS Toolbox
% Available Metrics: 'EUCLID', 'SEUCLID', 'CITYBLOCK', 'MAHAL', 'MINSKOWSKY'
% Available Methods for LINKAGE: 'SINGLE', 'COMPLETE', 'AVERAGE', 'CENTROID', 'WARD'
% Internal functions: CLONE, SUPPRESS, VER_EQ, EXTRACT
%
% vD	-> vector of distances (triu(D))
% Z	-> hierarchical clustering info using the k-NN algorithm
% I	-> inconsistent values (based on a depth meausre) of a cluster tree
% H	-> vector of line handles
% Cn	-> vector of cluster number for each cell
%

depth = 2;

vD = pdist(M,'euclid');
Z = linkage(vD,'centroid');
% I = inconsistent(Z,depth);
[H,Cn] = dendrogram(Z,0); % DO NOT work for 'WARD' linkage
ylabel('Inter-cell Affinities'); title('Network Dendrogram');


% SECONDARY INTERNAL FUNCTIONS %

function [h,T] = dendrogram(Z,p)
%DENDROGRAM Generate dendragram plot.
%   DENDROGRAM(Z) generates a dendrogram from the output matrix of
%   LINKAGE.  Z is a (M-1) by 3 matrix. M is the number of
%   observations in the original data. 
%
%   A dendrogram consists of many upsidedown U shape lines connecting
%   nodes in a hierachichale tree. Except for the WARD linkage (see
%   LINKAGE), the height of each U is the distance between the two
%   clusters to be connected at that time.  
%
%   DENDROGRAM(Z,P) generate a dendrogram with only the top P nodes.
%   When there are more than 30 nodes in the original data, the
%   dendrogram may look crowded. The default value of P is 30. If P =
%   0, then, every node will be displayed.
%
%   H = DENDROGRAM(...) returns a vector of line handles.
%
%   [H, T] = DENDROGRAM(...) also returns T, a vector of size M that
%   contains cluster number for each observation in the original data.
%
%   When bottom leaves are cutoff, some information are lost. T
%   supplies this lost information. For example, to find out which
%   observations are contained in node k of the dendrogram, use
%   FIND(T==k). 
%
%   When there are less than P observations in the original data, T
%   is the identical map, i.e. T = (1:M)'. Each node only contains
%   itself.
%
%   See also LINKAGE, PDIST, CLUSTER, CLUSTERDATA, INCONSISTENT.

%   ZP You, 3-10-98
%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 1.2 $

m = size(Z,1)+1;
if nargin < 2
   p = 30;
end

Z = transz(Z); % convert from m+k indexing to min(i,j) indexing.
T = (1:m)';

% if there are more than p node, dendrogram looks crowded, the following code
% will make the last p link nodes as the leaf node.
if (m > p) & (p ~= 0)
   
   Y = Z((m-p+1):end,:);
   
   R = Y(:,1:2);
   R = unique(R(:));
   Rlp = R(R<=p);
   Rgp = R(R>p);
   W(Rlp) = Rlp;
   W(Rgp) = setdiff(1:p, Rlp);
   W = W';
   T(R) = W(R);
   
   % computer all the leaf that each node (in the last 30 row) has 
   for i = 1:p
      c = R(i);
      T = clusternum(Z,T,W(c),c,m-p+1,0); % assign to it's leaves.
   end
   
   
   Y(:,1) = W(Y(:,1));
   Y(:,2) = W(Y(:,2));
   Z = Y;
   Z(:,3) = Z(:,3)-min(Z(:,3))*0.8; % this is to make the graph look more uniform.
   
   m = p; % reset the number of node to be 30 (row number = 29).
end

A = zeros(4,m-1);
B = A;
n = m;
X = 1:n;
Y = zeros(n,1);
r = Y;

% arrange Z into W so that there will be no crossing in the dendrogram.
W = zeros(size(Z));
W(1,:) = Z(1,:);

nsw = zeros(n,1); rsw = nsw;
nsw(Z(1,1:2)) = 1; rsw(1) = 1;
k = 2; s = 2;

while (k < n)
   i = s;
   while rsw(i) | ~any(nsw(Z(i,1:2)))
      if rsw(i) & i == s, s = s+1; end
      i = i+1;
   end
   
   W(k,:) = Z(i,:);
   nsw(Z(i,1:2)) = 1;
   rsw(i) = 1;
   if s == i, s = s+1; end
   k = k+1;
end

g = 1;
for k = 1:m-1 % initialize X
   i = W(k,1);
   if ~r(i), 
      X(i) = g; 
      g = g+1; 
      r(i) = 1;   
   end
   i = W(k,2);
   if ~r(i), 
      X(i) = g; 
      g = g+1; 
      r(i) = 1;   
   end
end
[u,v]=sort(X);
label = num2str(v');

for n=1:(m-1)
   i = Z(n,1); j = Z(n,2); w = Z(n,3);
   A(:,n) = [X(i) X(i) X(j) X(j)]';
   B(:,n) = [Y(i) w w Y(j)]';
   X(i) = (X(i)+X(j))/2; Y(i) = w;
end

% figure
set(gcf,'Position', [50, 50, 800, 500]);
h = plot(A,B,'b');
axis([0 m+2 0 max(Z(:,3))*1.05])
set(gca,'XTickLabel',[],'XTick',[],'box','off','xcolor','w');

if m <= 50,
   text((1:m)-0.2,zeros(m,1)-0.05*max(Z(:,3)),label);
end;

function T = clusternum(X, T, c, k, m, d)
% assign leaves under cluster c to c.

d = d+1;
n = m; flag = 0;
while n > 1
  n = n-1;
  if X(n,1) == k % node k is not a leave, it has subtrees 
     T = clusternum(X, T, c, k, n,d); % trace back left subtree 
     T = clusternum(X, T, c, X(n,2), n,d);
     flag = 1; break;
  end
end

n = size(X,1);
if flag == 0 & d ~= 1 % row m is leaf node.
   T(X(m,1)) = c;
   T(X(m,2)) = c;
end


function Y = pdist(X,s,t)
%PDIST Pairwise distance between observations.
%   Y = PDIST(X,METRIC) returns a vector which contains all the
%   distances between each pair of observations in X computed using
%   the given METRIC.  X is a M by N matrix, treated as M observations
%   of N variables. Since there are M*(M-1)/2 pairs of observations in
%   X, the size of Y is M*(M-1)/2 by 1.  The default metric is
%   'EUCLID'.  The available metrics are:
%
%      'euclid'    --- Euclidean metric
%      'seuclid'   --- Standardized Euclid metric
%      'cityblock' --- City Block metric
%      'mahal'     --- Mahalanobis metric
%      'minkowski' --- Minkowski metric
%
%   Y = PDIST(X, 'minkowski', p) specifies the exponents in the
%   'Minkowski' computation. When p is not given, p is set to 2.
%
%   The output Y is arranged in the order of ((1,2),(1,3),..., (1,M),
%   (2,3),...(2,M),.....(M-1,M)).  i.e. the upper right triangle of
%   the M by M square matrix. To get the distance between observation
%   i and observation j, either use the formula Y((i-1)*(M-i/2)+j-i)
%   or use the helper function Z = SQUAREFORM(Y), which will return a
%   M by M symmetric square matrix, with Z(i,j) equaling the distance
%   between observation i and observation j.
%
%   See also SQUAREFORM, LINKAGE

%   ZP You, 3-10-98
%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 1.2 $

if nargin >= 2
   if length(s) < 2
      error('Unrecognized metric');
   else 
      s = lower(s(1:2));
   end
else
   s = 'eu';
end

if s == 'mi' % Minkowski distance need a third argument
   if nargin < 3
      t = 2; 
   elseif t <= 0
      error('The third argument has to be positive.');
   end
end

[m, n] = size(X);

if m < 2
   error('The first argument has to be a numerical matrix with at least two rows');
end

p = (m-1):-1:2;
I = zeros(m*(m-1)/2,1);
I(cumsum([1 p])) = 1;
I = cumsum(I);
J = ones(m*(m-1)/2,1);
J(cumsum(p)+1) = 2-p;
J(1)=2;
J = cumsum(J);

Y = (X(I,:)-X(J,:))';
I = []; J = []; p = [];  % no need for I J p any more.

switch s
case 'eu' % Euclidean
   Y = sum(Y.^2);
   Y = sqrt(Y);
case 'se' % Standadized Euclidean
   D = diag(var(X));
   Y = sum(D*(Y.^2));
   Y = sqrt(Y);
case 'ci' % City Block
   Y = sum(abs(Y));
case 'ma' % Mahalanobis
   v = inv(cov(X));
   Y = sqrt(sum((v*Y).*Y));
case 'mi' % Minkowski
   Y = sum(abs(Y).^t).^(1/t);
otherwise
   error('no such method.');
end


function Z = linkage(Y, method)
%LINKAGE Hierarchical cluster information.
%   LINKAGE(Y) computes the hierarchical cluster information, using the
%   single linkage algorithm, from a given distance matrix Y generated
%   by PDIST. Y is also commonly known as similarity or
%   dissimilarity matrix.
%
%   LINKAGE(Y, method) computes the hierarchical cluster information using
%   the specified algorithm. The available methods are:
%
%      'single'   --- nearest distance
%      'complete' --- furthest distance
%      'average'  --- average distance
%      'centroid' --- center of mass distance
%      'ward'     --- inner squared distance
%
%   Cluster information will be returned in the matrix Z with size m-1
%   by 3.  Column 1 and 2 of Z contain cluster indices linked in pairs
%   to form a binary tree. The leaf nodes are numbered from 1 to
%   m. They are the singleton clusters from which all higher clusters
%   are built. Each newly-formed cluster, corresponding to Z(i,:), is
%   assigned the index m+i, where m is the total number of initial
%   leaves. Z(i,1:2) contains the indices of the two component
%   clusters which form cluster m+i. There are n-1 higher clusters
%   which correspond to the interior nodes of the output clustering
%   tree. Z(i,3) contains the corresponding linkage distances between
%   the two clusters which are merged in Z(i,:), e.g. if there are
%   total of 30 initial nodes, and at step 12, cluster 5 and cluster 7
%   are combined and their distance at this time is 1.5, then row 12
%   of Z will be (5,7,1.5). The newly formed cluster will have an
%   index 12+30=42. If cluster 42 shows up in a latter row, that means
%   this newly formed cluster is being combined again into some bigger
%   cluster.
%
%   See also PDIST, INCONSISTENT, COPHENET, DENDROGRAM, CLUSTER, CLUSTERDATA

%   ZP You, 3-10-98
%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 1.4 $

[k, n] = size(Y);

if n < 3
  error('You have to have at least three distances to do a linkage.');
end
  

m = (1+sqrt(1+8*n))/2;
if k ~= 1 | m ~= fix(m)
  error('The first input has to match the output of the PDIST function in size.');   
end

if nargin == 1 % set default switch to be 's' 
   method = 'si';
end

if length(method) < 2
   error('The switch given by the second argument is not defined.');
end

method = lower(method(1:2)); % simplify the switch string.

Z = zeros(m-1,3); % allocate the output matrix.

% during updating clusters, cluster index is constantly changing, R is
% a index vector mapping the original index to the current (row
% column) index in X.  N denotes how many points are contained in each
% cluster.

N = zeros(1,2*m-1);
N(1:m) = 1;
n = m; % since m is changing, we need to save m in n. 
R = 1:n;

if method == 'ce'  % square the X so that it is easier to update.
   Y = Y .* Y;
elseif method == 'wa'
   Y = Y .* Y/2;
end

for s = 1:(n-1)
   if method == 'av'
      p = (m-1):-1:2;
      I = zeros(m*(m-1)/2,1);
      I(cumsum([1 p])) = 1;
      I = cumsum(I);
      J = ones(m*(m-1)/2,1);
      J(cumsum(p)+1) = 2-p;
      J(1)=2;
      J = cumsum(J);
      W = N(R(I)).*N(R(J));
      X = Y./W;   
   else
      X = Y;
   end
   
   [v, k] = min(X);
   if method == 'ce'
      v = sqrt(v);
   end
   
   i = floor(m+1/2-sqrt(m^2-m+1/4-2*(k-1)));
   j = k - (i-1)*(m-i/2)+i;
   
   Z(s,:) = [R(i) R(j) v]; % update one more row to the output matrix A
   
   % update X, in order to vectorize the computation, we need to compute
   % all the index corresponds to cluster i and j in X, denoted by I and J.
   I1 = 1:(i-1); I2 = (i+1):(j-1); I3 = (j+1):m; % these are temp variables.
   U = [I1 I2 I3];
   I = [I1.*(m-(I1+1)/2)-m+i i*(m-(i+1)/2)-m+I2 i*(m-(i+1)/2)-m+I3];
   J = [I1.*(m-(I1+1)/2)-m+j I2.*(m-(I2+1)/2)-m+j j*(m-(j+1)/2)-m+I3];
   
   switch method
   case 'si' %single linkage
      Y(I) = min(Y(I),Y(J));
   case 'av' % average linkage
      Y(I) = Y(I) + Y(J);
   case 'co' %complete linkage
      Y(I) = max(Y(I),Y(J));
   case 'ce' % centroid linkage
      K = N(R(i))+N(R(j));
      Y(I) = (N(R(i)).*Y(I)+N(R(j)).*Y(J)-(N(R(i)).*N(R(j))*v^2)./K)./K;
   case 'wa'
      Y(I) = ((N(R(U))+N(R(i))).*Y(I) + (N(R(U))+N(R(j))).*Y(J) - ...
	  N(R(U))*v)./(N(R(i))+N(R(j))+N(R(U)));
   otherwise error('method not recognized.');
   end
   J = [J i*(m-(i+1)/2)-m+j];
   Y(J) = []; % no need for the cluster information about j.
   
   % update m, N, R
   m = m-1; 
   N(n+s) = N(R(i)) + N(R(j));
   R(i) = n+s;
   R(j:(n-1))=R((j+1):n); 
end

if method == 'wa'
   Z(:,3) = sqrt(Z(:,3));
end

function Z = transz(Z)
%TRANSZ Translate output of LINKAGE into another format.
%   This is a helper function used by DENDROGRAM and COPHENET.

%   In LINKAGE, when a new cluster is formed from cluster i & j, it is
%   easier for the latter computation to name the newly formed cluster
%   min(i,j). However, this definition makes it hard to understand
%   the linkage information. We choose to give the newly formed
%   cluster a cluster index M+k, where M is the number of original
%   observation, and k means that this new cluster is the kth cluster
%   to be formmed. This helper function converts the M+k indexing into
%   min(i,j) indexing.

%   ZP You, 3-10-98
%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 1.1 $

m = size(Z,1)+1;

for i = 1:(m-1)
   if Z(i,1) > m      
      Z(i,1) = traceback(Z,Z(i,1));
   end
   if Z(i,2) > m
      Z(i,2) = traceback(Z,Z(i,2));
   end
   if Z(i,1) > Z(i,2)
      Z(i,1:2) = Z(i,[2 1]);
   end
end

function a = traceback(Z,b)

m = size(Z,1)+1;

if Z(b-m,1) > m
   a = traceback(Z,Z(b-m,1));
else
   a = Z(b-m,1);
end
if Z(b-m,2) > m
   c = traceback(Z,Z(b-m,2)); 
   else
   c = Z(b-m,2);
end

a = min(a,c);
              