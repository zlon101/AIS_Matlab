function [E,Nc,mCe,mC,T] = anal1(M,D,st,s);

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

% disp(sprintf('FUNCTION: ANALYSIS'));
% disp(sprintf('Copyright by Leandro de Castro, March - 2000\n'));
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
if N > 5,
   figure(3); clf; bar(E(:,3)); title('Number of Clusters (Valleys)');
else,
   mC = 1:1:5; Nc = 5; mCe = M; T = mC;
   break;
end;

% Determining the Number and Separation of clusters by using the E matrix
figure(4); clf; draw_net(M); hold on;
T = clusterdata(M,5); Nc = max(T);
mC = zeros(Nc,length(T)); mCe = []; max_aux = 0;
for i=1:Nc
   aux = find(T==i)'; comp = length(aux);
   mC(i,1:comp) = aux;
   if comp > max_aux,
      max_aux = comp; 
   end;
   if comp > 1,
      mCe(i,:) = mean(M(aux,:));
   else,
      mCe(i,:) = M(aux,:);
   end;
   d_net(M(aux,:),1,2,1);
end;
mC = mC(:,1:max_aux); title('Final Network Structure');
plot(mCe(:,1),mCe(:,2),'r*'); hold off; 
% Then separate cells by clusters, determine centroid and plot network

% Call function Dendro to draw the Dendrogram
figure(5); [Z,H,Cn] = dendro(M);

% ------------------------------------- %
%           End of Main Function        %
% ------------------------------------- %

% SECONDARY INTERNAL FUNCTIONS %

% Function Draw Network
function [D,I] = draw_net(M,Da,Na,td);
if nargin == 1,
   Na = 2; td = 0.5; Da = 0;
elseif nargin == 2
   Na = 2; td = 1.5;
end;
[N1,L] = size(M);
if L < 2, disp('Improper number of columns'); break; end;
if Na > N1-1, disp('Improper number of arcs'); Na = 1; end;
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


% Function Draw Network
function [D,I] = d_net(M,Da,Na,td);
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
   % a = num2str(i); text(M(i,1)+val,M(i,2)+val,a); % Ab indexes
   for j=1:Na,
      if D(i,I(j,i)) < td & Da == 1,
         line([M(i,1),M(I(j,i),1)],[M(i,2),M(I(j,i),2)]);
      end;
   end;
end;
% End Function DRAW_NET


