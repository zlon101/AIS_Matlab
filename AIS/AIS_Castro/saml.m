function [ab,vE,E] = saml(N);
%
% Ph.D. Thesis
% Leandro Nunes de Castro
% October, 1999.
% Function generates n-bit binary strings antibody most uniformly 
% distributed over the search space
% Simulated Annealing (SA)
% Characteristics:
% * The fitness of each individual is a function of the others (cooperation)
% * Mutation of the whole population
%
% Auxiliar functions: CADEIA, MATCH, HYPERMUT
%
% [ab,vE,Tempo,Flops] = sa1(ab,niter);
% ab		-> population of antibodies
% vE		-> annealing schedule
% niter  -> number of iterations
%
% Definitions:
% D   -> percentile amount of non-repeated individuals
% F   -> percentile HD of the population
% E   -> Energy measure (percentile)
%

if nargin == 1,
   [ab] = ones(N,22);
   niter = 500 * N;
end;
figure(1); clf; plotfun(ab);
SC = 100;
disp(sprintf('Maximal number of iteration steps: %d',niter));
disp(sprintf('Stopping criterion: %.1f',SC));
disp(sprintf('Size of vectors: [%d,%d]',N,22));
disp(sprintf('Press any key to continue...'));pause;

vE = []; it = 1; [N,s1] = size(ab);
% Simulated Annealing parameters
alfa = 0.8;        % Geometrical cooling
T = 1;             % Initial temperature

f0 = flops; t0 = clock;
% First Step
for i = 1:N,
   [ms(i,:),d(i)] = match(ab,ab(i,:));
end;
f = (100/(N*s1)) .* sum(ms);
D = (100/N) .* sum(d == (N-1));
F = 2/N .* sum(f); E = (F + D)/2;
Ea = E; lr = 0.1;

% Iterative Process
while it <= niter & E ~= SC,
   % Simulated Annealing
   aux = hypermut(ab,lr);
   for i = 1:N,
      [ms(i,:),d(i)] = match(aux,aux(i,:));
   end;
   f = (100/(N*s1)) .* sum(ms); F = 2/N .* sum(f);
   D = (100/N) .* sum(d == (N-1));
   E = (F + D)/2; dE = Ea - E;
   if dE < 0,
      ab = aux;
   elseif dE > 0,
      if rand < exp(-dE/T), % Probabilistic behavior
         %disp('** Annealing **');
         ab = aux;
      else,
         E = Ea;
      end;
   elseif rem(it,5) == 0,   % Identify Steady-State
      % disp('** Cooling **');
      T = alfa * T;
   end;
   if rem(it,20) == 0,
      disp(sprintf('It: %d	T: %2.3f	LR: %2.4f	E: %2.2f',it,T,lr,E));
      %if lr > 5e-3,
      %   lr = 0.9*lr;
      %end;
   end;
   if rem(it,50) == 0,
      figure(2); clf; plotfun(ab); drawnow;
   end;
   it = it + 1; Ea = E; vE = [vE;E];
end;
disp(sprintf('It: %d	T: %2.3f	LR: %2.4f	E: %2.2f',it,T,lr,E));


% clf; plot(vE); title('Annealing Evolution'); ylabel('Energy'); xlabel('Iterations');
% axis([0 length(vE)+0.02*length(vE) min(vE)-0.02*min(vE) vE(end)+0.02*min(vE)]);
 
% --------------------- %
% INTERNAL SUBFUNCTIONS %
% --------------------- %
 
% Function DECODE
function x = decode(v);
if min(min(v)) == -1, v = hardlim(v); end;
v = fliplr(v); s = size(v);
aux = 0:1:21; aux = ones(s(1),1)*aux;
x1 = sum((v.*2.^aux)');
x = -1 + x1 .* (2 / 4194303);
% End Function DECODE

% Function CADEIA
function [ab,ag] = cadeia(n1,s1,n2,s2,bip)
if nargin == 2,
   n2 = n1; s2 = s1; bip = 1;
elseif nargin == 4,
   bip = 1;
end;
% Antibody (Ab) chains
ab = 2 .* rand(n1,s1) - 1;
% ab = randn(n1,s1);
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
% End function CADEIA

% Function MATCH
function [ms,D] = match(ab,ag,comp)
if nargin == 2,
   comp = 0;                     % Hamming distance
end;
msc = [];   % ms complement
% Converting bipolar (-1,+1) strings to binary ones (0,+1)
ab = hardlim(ab);
ag = hardlim(ag);
% Using the XOR operator for calculating the match score
[n1,s1] = size(ab);
% Batch Mode
ag = ones(n1,1) * ag;					% Multiply the Antigen
mXOR = xor(ab,ag);
ms = sum(mXOR');
D = sum(ms > 0);
msc = 1 - ms;
if comp == 1,
   ms = msc;
end;
% End function MATCH

% Function HYPERMUT
function [v] = hypermut(v,pm),
[npop,s1] = size(v);
mat = rand(npop,s1) <= pm;
v = v - 2 .* (v .* mat);
% End function HYPERMUT

% Function PLOTFUN
function [] = plotfun(ab);
x = -1:.01:1;
f = 'x .* sin(10 * pi .* x) + 1'; y = eval(f);
plot(x,y,'k'); hold on;
x = decode(ab);
f = 'x .* sin(10 * pi .* x) + 1'; y = eval(f);
plot(x,y,'ro'); hold off;
% End function PLOTFUN
