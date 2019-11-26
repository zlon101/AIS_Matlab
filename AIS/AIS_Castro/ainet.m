function [M,D] = inet(Ag,ts,n,N,gen,qi,tp)

%
% Ph.D. Thesis
% Leandro Nunes de Castro
% February, 2000
% Immune Network (iNet) - Description in iNet.doc
% Data normalization over [0,1] required 
% Obs.: for simplicity of comprehension we chose non-complementary vectors
% Further extension: complementary vectors
% Secondary Functions: RUN_INET, DIST, CENTROID
% Internal functions: CLONE, SUPPRESS, VER_EQ, EXTRACT, PLOTVET1, DRAW_NET, NORMA
%
% function [M,D] = inet(Ag,Ab,Agt,n,N,gen)
% M     -> memory cells matrix
% D     -> distance matrix for M
% Ag    -> antigens (training patterns)
% Ab    -> network antibodies
% n     -> no. of best-matching cells taken for each Ag (Selection)
% N     -> clone number multiplier
% gen   -> maximum number of generations
%
% L     -> Ag and Ab length
% N1    -> no. of antibodies (constructive)
% N2    -> no. of antigens
% Nc    -> no. of clones to be generated
% D     -> Ag-Ab affinity vector
% Do    -> D sorted in ascending order
% mi    -> learning (hypermutation) rate (default: 4.0)
% qi    -> percentile amount of clones to be Re-selected
% ts    -> suppression threshold (default: 0.001)
% tp    -> pruning threshold
% D1    -> idiotypic affinity matrix [N1,N1]
% vbD   -> vector best affinity for each Ag
% nR    -> no. of Ab to be re-selected
%
% Adequate parameters for solving the SPIRAL task: td = 0.1; ts = 0.04; gen = 100;

[N2,L] = size(Ag);
N1 = 10; Ab = rand(N1,L); mi = 4.0; Agt = rand(N1,L); sc = 0.01;
   
% disp(sprintf('Suppression threshold: %.2f',ts));
disp(sprintf('Pruning threshold: %.2f',tp/100));
disp(sprintf('Number of best matching cells to be selected: %d',n));
disp(sprintf('Percentile amount of clones to be re-selected: %d%%',100*qi));
disp(sprintf('Number of generations: %d',gen));
disp(sprintf('Population (to be recognized) size: [%d,%d]',N2,L));
disp('Press any key to continue...'); pause;

it = 0; avD = 1000;
RUNNET = 0;
while it < gen & avD > sc,
   vbD = []; M = []; i = 1;
   vet = randperm(N2);
   while i <= N2,
      % Ag-Ab Affinity (Match-Function: Euclidian distance)
      D = dist(Ab,Ag(vet(i),:)');
      [Dn,I] = sort(norma(D));
      Nc = floor(N-Dn(1:n,:)*N);
      
      % Clone & Affinity Maturation
      [C,Cag,Cmi] = clone(Ab,Ag(vet(i),:),mi,D,I,Nc);
      C = C - Cmi.*(C-Cag);
      
      % Re-Selection
      D = dist(C,Ag(vet(i),:)');
      [Dn,I] = sort(D);
      nR = round(qi*size(C,1));
      m = C(I(1:nR),:); % 1 clone for each Ag
      D = D(I(1:nR));   % new affinities
      
      % Network pruning (Natural Death)
      Ip = find(D > tp);
      m = extract(m,Ip); D = extract(D,Ip);
      
      % Suppression (Idiotypic Network)
      [m,D1] = suppress(m,ts);
      
      % General parameters
      minD = min(D); [vbD] = [vbD; minD];
      Cs = size(m,1);
      M = [M; m];  % memory matrix
      % disp(sprintf('Ag: %d	minD: %.4f	Cs: %d	mi: %.4f',i,minD,Cs,mi));
      
      % Utilizar uma função p/ controlar o valor de mi
      %mi = 0.9 * mi; % decreasing mi
      i = i + 1;
   end;
   % Search for similaritites among clusters
   [M,D] = suppress(M,ts); 
   
   % Re-build Ab repertoire
   Ab = [M;rand(N1,L)];
   avD = [avD;mean(vbD)];
   it = it + 1;
   
   disp(sprintf('It: %d	avD: %f	Net size: %d',it,avD(end),size(M,1)));
end;

% Clustering
%if RUNNET == 1,
%   Itrain  = run_inet(M,Ag);
%   Itest = run_inet(M,Agt);
%end;

% Drawing Network & Plotting Results
% figure(2); plot(avD(2:end)); title('Average Affinities'); xlabel('Generations');

% SECONDARY FUNCTIONS %

function [C,Cag,Cmi] = clone(Ab,ag,mi,D,I,Nc);
% C   -> clones (from greater to smaller affinity)
% Cag -> clones of Ag
% Cmi -> clones of mi
% S   -> selected antibodies
% Obs.: Cag and Cmi are necessary for the updating procedure
%       The original cell is mantained
[N1,L] = size(Ab); [n,N2] = size(Nc);
% C = []; Cmi = []; Cag = [];
C = Ab(I(1),:); Cmi = ones(1,L); Cag = C; % Maintenance of the fittest cell before maturation
for i=1:n,
   vones = ones(Nc(i),1);
   C = [C; vones * Ab(I(i),:)];
   Cag = [Cag; vones * ag];
   Cmi = [Cmi; rand(Nc(i),L) .* D(I(i)) .* mi];
end;

% Function suppress self-recognizing and non-stimulated Ab from Memory (M)
function [M,D1] = suppress(M,ts);
% M   -> memory matrix
% D1  -> idiotypic affinity matrix
D1 = dist(M,M');
aux = triu(D1,1);
[Is,Js] = find(aux>0 & aux<ts);
if ~isempty(Is),
   Is = ver_eq(Is);
   M = extract(M,Is);
   % D1 = extract(D1,Is);
end;
D1 = dist(M,M');

% Search for repeated indexes
function [Is] = ver_eq(I);
l = length(I); Is = [];
if l > 1,
   for i=1:l-1,
      aux = I(i);
      auxI = I(i+1:end);
      el = find(auxI == aux);
      if isempty(el),
         Is = [Is,aux];
      end;
   end;
   Is = [Is,I(end)];
else,
   Is = I;
end;

% Function Extracts lines from M indexed by I
function [M] = extract(M,I);
Maux = zeros(size(M));
Maux(I,:) = M(I,:);
M = M - Maux;
[I] = find(M(:,1)~=0);
M = M(I,:);

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
