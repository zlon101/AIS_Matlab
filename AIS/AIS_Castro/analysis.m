function [E,bE,Nc,mCe,mC,T,U] = analysis(M,D,st,s);

%
% Ph.D. Thesis
% Copyright by Leandro Nunes de Castro
% March, 2000
% Immune Network (iNet) - Description in iNet.doc
% Function determines the Minimal Spanning Tree (MST) of the aiNet
% Number and Members of each Cluster
% Source: Discrete Mathematics, R. Johnsonbaugh, pp.401 (1997)
% Secondary Functions: DRAW_NET, UNIQUE, DENDRO, SEPARATE
%
% function [E,bE,Nc,mCe,mC,T,U] = analysis(M,D,s,st);
% E    -> set of edges and the distance between them: [e1,e2,d(e1,e2)]
% bE   -> edges that separate clusters (Col 3 is the index of the cutting line of E)
% Nc   -> number of clusters
% mCe  -> matrix containing the centroids of each cluster
% mC   -> matrix containing the cells for each cluster
% T    -> matrix labeling each cell to a cluster i: [1xi] 
% U    -> fuzzy membership matrix: [cellxcluster]
% M    -> matrix of memory cell coordinates
% D    -> distance matrix among cells
% s    -> starting edge
% st   -> std threshold
%
% v(i) = 1 if vertex i has been added to MST
% v(i) = 0 if vertex i has not been added to MST
% bE   -> breaking edge
% mC   -> centroid matrix
%

disp(sprintf('FUNCTION: ANALYSIS'));
disp(sprintf('Copyright by Leandro de Castro, March - 2000\n'));
disp(sprintf('** This Function Determines **:	\n1. Minimal Spanning Tree\n2. Cluster Analysis\n3. Centroid and Cluster Membership\n4. Dendrogram'));
if nargin == 2, 
   s = 1; st = 2;
elseif nargin == 3,
   s = 1;
end;
N = length(D);
if s > N, disp('Improper initial vertex'); s = 1; end;

% Drawing the Minimal Spanning Tree (MST)
v = [zeros(N,1)]; v(s) = 1; E = []; 
figure(2); clf; hold on; draw_net(M);
for i=1:N-1,
   menor = 1e3;
   for j=1:N,
      if v(j) == 1,  % j is a vertex in MST
         for k=1:N,
            if v(k) == 0 & D(j,k) < menor,
               add_vertex = k;
               e = [j,k];
               menor = D(j,k);
            end;
         end;
      end;
   end;
   v(add_vertex) = 1;
   E = [E;e,D(e(1),e(2))];
   line([M(E(i,1),1),M(E(i,2),1)],[M(E(i,1),2),M(E(i,2),2)]);
end;
title('Minimal Spanning Tree'); hold off;

% Visualizing the Cluster Separation based on matrix E
%T = dist(M(E(end,2),:),M(E(1,1),:)'); T = [T;E(:,3)];
%figure(2); clf; bar(T); title('Number of Clusters (Peaks)');
figure(3); clf; bar(E(:,3)); title('Number of Clusters (Peaks + 1)');

% Determining the Number and Separation of clusters by using the E matrix
figure(4); clf; draw_net(M); hold on;
[bE,Nc,mCe,mC,T] = separate(E,N,M,st);
U = dist(M,mCe'); % U = normal(U); 
% Correct U: 18/03/00, try to apply a membership function

% Test CLUSTERDATA
% T = clusterdata(M,0.8);
% Then separate cells by clusters, determine centroid and plot network

% Call function Dendro to draw the Dendrogram
figure(5); [Z,H,Cn] = dendro(M);

% ------------------------------------- %
%           End of Main Function        %
% ------------------------------------- %

% SECONDARY INTERNAL FUNCTIONS %

% Determining the Number and Separation of clusters by using the E matrix
function [bE,Nc,mCe,mC,T] = separate(E,N,M,st);
vE = []; vC = []; Nc = 1; bE = [];
flag = 0; i = 1;
while i < N-1,
   vE = [vE;E(i,3)]; flag = 0;
   mE = mean(vE); stdE = std(vE);
   % Ratio between d(x,y) and the respective averages
   if E(i+1,3) > (mE + st*stdE),
      vE = [];
      Nc = Nc + 1; i = i + 1;
      bE = [bE;E(i,1:2),i];
   end;
   i = i + 1;
end;
if Nc == 1,
   disp(sprintf('Attention: Single Cluster\nReduce the STD multiplier'));
   disp(sprintf('Actual value: %f',st));
   mCe = []; mC = []; T = ones(1,N);
   break;
end;

% Centroids (mCe) and Members of each cluster (mC)
ini = 1; mCe = []; mC = zeros(Nc,N);  maxim = 0; T = zeros(1,N);
% Draw Network and Centroid, Cluster by Cluster
for i=1:Nc-1,
   if ini < bE(i,3)-1,
      Nodes = reshape(E(ini:(bE(i,3)-1),1:2),1,(bE(i,3)-ini)*2);
   else,
      Nodes = E(ini,1:2);
   end;
   Nodes = unique(Nodes);
   T(Nodes) = i;
   lN = length(Nodes);
   mC(i,1:lN) = Nodes;
   maxim = max(maxim,lN);
   ini = bE(i,3) + 1;
   mCe = [mCe;mean(M(Nodes,:))];
   draw(M(Nodes,:),1,2,2); plot(mCe(i,1),mCe(i,2),'r*');
end;
Nr = size(E,1) + 1;
if Nr > ini,  % Last rows (cluster) of mC and mCe
   Nodes = reshape(E(ini:Nr,1:2),1,(Nr-ini+1)*2);
   Nodes = unique(Nodes);
   mC(end,1:length(Nodes)) = sort(Nodes);
   maxim = max(maxim,length(Nodes));
   mCe = [mCe;mean(M(Nodes,:))];
else,         % last row (cluster) of mC and mCe
   if T(E(end,1)) == 0 & T(E(end,2)) == 0,
      Nodes = E(end,1:2); T(E(end,1:2)) = Nc;
      mCe = [mCe;mean(M(Nodes,:))];
   elseif T(E(end,1)) == 0 & T(E(end,2)) ~= 0,
      Nodes = E(end,1); T(E(end,1)) = Nc;
      mCe = [mCe;M(Nodes,:)];
   elseif T(E(end,1)) ~= 0 & T(E(end,2)) == 0,
      Nodes = E(end,2); T(E(end,2)) = Nc;
      mCe = [mCe;M(Nodes,:)];
   end;
   mC(end,1:2) = sort(Nodes);
end;
draw(M(Nodes,:),1,2,2); plot(mCe(end,1),mCe(end,2),'r*');
mC = mC(:,1:maxim); hold off;

% Function draw the network
function [D,I] = draw(M,Da,Na,td);
if nargin == 1,
   Na = 2; td = 0.5; Da = 0;
elseif nargin == 2,
   Na = 2; td = 1.5;
end;
% Only draw the connection to the closest cell and to the cluster centroid
[N1,L] = size(M);
if L < 2, disp('Improper number of columns'); break; end;
if Na > N1-1, Na = N1-1; end;
axis([-0.1 1.1 -0.1 1.1]);

% Determines the affinity between each Ab
D = dist(M,M');
[aux,I] = sort(D);
I = I(2:N1,:); % Eliminate itself
% val = 0.01*max(max(M));
for i=1:N1,
   % a = num2str(i); text(M(i,1)+val,M(i,2)+val,a); % Ab indexes
   for j=1:Na,
      if D(i,I(j,i)) < td & Da == 1,
         line([M(i,1),M(I(j,i),1)],[M(i,2),M(I(j,i),2)]);
      end;
   end;
   % line([M(i,1),mCe
end;


% Function Draw Network
function [D,I] = draw_net(M,Da,Na,td);
if nargin == 1,
   Na = 2; td = 0.5; Da = 0;
elseif nargin == 2
   Na = 2; td = 1.5;
end;
[N1,L] = size(M);
if L < 2, disp('Improper number of columns'); break; end;
if Na > N1-1, disp('Improper number of arcs'); Na = N1-1; end;
% For visualization purposes it is not necessary to identify all coordinates
plot(M(:,1),M(:,2),'ko'); drawnow; hold on;
axis([-0.1 1.1 -0.1 1.1]);
% Determines the affinity between each Ab
D = dist(M,M');
[aux,I] = sort(D);
I = I(2:N1,:); % Eliminate itself
val = 0.01*max(max(M));
for i=1:N1,
   a = num2str(i); text(M(i,1)+val,M(i,2)+val,a); % Ab indexes
   for j=1:Na,
      if D(i,I(j,i)) < td & Da == 1,
         line([M(i,1),M(I(j,i),1)],[M(i,2),M(I(j,i),2)]);
      end;
   end;
end;
% End Function DRAW_NET

