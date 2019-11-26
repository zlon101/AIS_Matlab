function [M] = imalgchar(P,X,gen,n,pm,per);

%
% Ph.D. Thesis
% Leandro Nunes de Castro
% November, 1999.
% Immune Algorithm - New evolutionary strategy inspired in the Immune System
% Operations: Hypermutation, Editing, Selection
% Each clone has size proportional to its affinity
%
% function [M] = pattern(P,X,gen,n,pm,per);
% M      -> memory matrix
% P		-> population of size N x L
% X		-> patterns to be recognized np x L
% gen	   -> number of generations
% n		-> number of clones
% pm	   -> hypermutation probability
% per  	-> percentile of the population to suffer random reshuffle
% 
% T		-> temporary population
% M		-> Memory matrix (functionally disconnected)
% Tips: proportional clone sizes with very high sizes, e.g., fat = L
%

if nargin == 2,
   gen = 100; n = round(size(P,1)/2); pm = 0.05; per = 0.0; fat = 10;
   % gen = 200; n = size(P,1); pm = 0.1; per = 0.0; fat = 35;
end;
while n <= 0,
   n = input('n has to be at least one. Type a new value for n: ');
end;

[N,L] = size(P); it = 0;
np = size(X,1); PRINT = 1;
% Hypermutation controlling parameters
pma = pm; itpm = 10; pmr = 0.8;
mfit = []; vpm = []; menor = 1;
M = cadeia(np,L+1);
disp(sprintf('Population size: %d',N));
disp(sprintf('Memory matrix size: [%d,%d]',np,L));
disp(sprintf('Maximum number of generations: %d',gen));
imprime(PRINT,12,10,X,it,1,1); title('Pattern to be recognized');
imprime(PRINT,12,10,M(1,2:end),it,1,2); title('Initial memory matrix');
disp('Press any key to continue...'); pause;
M = hardlim(M); M(1,1) = L; % Transform into a binary matrix

% Generations
while it < gen & menor > 0,
   T = []; k = 0; vet = randperm(np); vfit = [];	 vind = [];		% Assincronous
   while k < np,
      k = k+1; fit = []; % i = vet(k); 
      [fit,mXOR] = match(P,X(k,:),0); 
      [v(k,:),ind(k,:)] = sort(fit);
      
      % Reproduction
      [T,pcs] = reprod(n,fat,N,ind(k,:),P,T);
      
      % Hypermutation
      Ta = rand(size(T,1),L) <= pm;
      T = T - 2 .* (T.*Ta) + Ta;    % 0,1 mutation
      T(pcs,:) = P(ind(k,1:n),:);	% keep the previous best individuals
      
      % Re-selection
      [fit,mXOR] = match(T,X(k,:),0);
      pcs = [0 pcs];
      for i=1:n,
         [out(i),bcs(i)] = min(fit(pcs(i)+1:pcs(i+1)));		% Minimization problem
         bcs(i) = bcs(i) + pcs(i);
      end;
      P(ind(k,1:n),:) = T(bcs,:);
      
      % Memory Assignment & Evaluation
      [b,indb] = min(fit);
      if b < M(k,1),
         M(k,1) = b; M(k,2:end) = T(indb,:);
         menor(k) = b;
      else,
         menor(k) = M(k,1);
      end;
      
      % Editing (Repertoire shift)
      nedit = round(per*N);
      if nedit > 0,
         P(ind(k,N-nedit-np+1:N-np),:) = cadeia(nedit,L,0,0,0);
      end;
   end;
   P(ind(:,1),:) = M(:,2:end);
   
   % Hypermutation control
   [pm] = pmcont(pm,pma,pmr,it,itpm);
   
   menor = sum(menor); vfit = [vfit,menor]; it = it + 1;
   disp(sprintf('It.: %d	pm: %.3f	F: %2.4f',it,pm,menor));
   Mem = hardlims(M - 0.1);
   imprime(PRINT,12,10,Mem(:,2:end),it,5,2);
end; % end while
imprime(PRINT,12,10,M(:,2:end),it,1,2);

% --------------------- %
% INTERNAL SUBFUNCTIONS
% --------------------- %

% Function plot pictures
function [] = imprime(PRINT,res_lin,res_col,P,it,mit,fn);
if PRINT == 1,
   if rem(it,mit) == 0,
      fig(res_lin,res_col,0,fn,P);
   end;
end;

function [mat] = fig(L,C,div,fn,X)
figure(fn); clf; hold on;
for i = 1:size(X,1),
   mat = reshape(X(i,:),C,L)';
   if (div == 1),
      subplot(2,4,i);
   end;
   image(mat*15);axis('square');axis('off');
end;
drawnow; hold off;

% Function match bipolar strings
function [ms,mXOR] = match(ab,ag,comp)
if nargin == 2,
   comp = 0;                     % Hamming distance
end;
msc = [];   % ms complement
if min(min(ag)) == -1,
   ag = hardlim(ag);
end;
% Using the XOR operator for calculating the match score
[n1,s1] = size(ab);
ag = ones(n1,1) * ag;					% Multiply the Antigen
mXOR = xor(ab,ag);
ms = sum(mXOR');
msc = 1 - ms;
if comp == 1,
   ms = msc;
end;

% Reproduction
function [T,pcs] = reprod(n,fat,N,ind,P,T);
% n		-> number of clones
% fat	-> multiplying factor
% ind	-> best individuals
% T		-> temporary population
% pcs	-> final position of each clone
if n == 1,
   pcs = N;
   T = ones(N,1) * P(ind(1),:);
else,
   for i=1:n,
      cs(i) = round(fat*N/i);
      % cs(i) = round(fat*N);
      pcs(i) = sum(cs);
      T = [T; ones(cs(i),1) * P(ind(i),:)];
   end;
end;

% Control of pm
function [pm] = pmcont(pm,pma,pmr,it,itpm);
% pma	-> initial value
% pmr	-> control rate
% itpm	-> iterations for restoring
if rem(it,itpm) == 0,
   pm = pm * pmr;
   if rem(it,10*itpm) == 0,
      pm = pma;
   end;
end;

function [ab,ag] = cadeia(n1,s1,n2,s2,bip)
if nargin == 2,
   n2 = n1; s2 = s1; bip = 1;
elseif nargin == 4,
   bip = 1;
end;
% Antibody (Ab) chains
ab = 2 .* rand(n1,s1) - 1;
if bip == 1,
   ab = hardlims(ab);
else,
   ab = hardlim(ab);
end;
% Antigen (Ag) chains
ag = 2 .* rand(n2,s2) - 1;
if bip == 1,
   ag = hardlims(ag);
else,
   ag = hardlim(ag);
end;