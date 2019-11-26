function [ab,vE,vRb,E,Rb] = sareal(ab,niter);
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
% Internal Subfunctions: NOISE, UNIT_NORM, ENER, PLOTVET
%
% [ab,vE,Tempo,Flops] = sareal(ab,niter);
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
   % N = 2; s = 2; [ab] = randn(N,s); 
   [N,s] = size(ab);
   ab = unit_norm(ab);
   niter = 2000 * N;
end;
figure(1); plotvet(ab); plot([0,0],[0,0],'ro'); drawnow; title('Initial Directions');
SC = 99.9;
disp(sprintf('Maximal number of iteration steps: %d',niter));
disp(sprintf('Stopping criterion: %.1f',SC));
disp(sprintf('Size of vectors: [%d,%d]',N,s));
disp(sprintf('Press any key to continue...'));pause;
vE = []; vRb = []; it = 0; [N,s1] = size(ab);
% Simulated Annealing parameters
alfa = 0.9;        % Geometrical cooling
T = 1;             % Initial temperature

% First Step
f0 = flops; t0 = clock;
[Ea,Rb] = ener(ab); lr = 1;

% Iterative Process
while it <= niter & Rb < SC,
   % Simulated Annealing
   aux = noise(ab,lr);
   [E,Rb] = ener(aux); dE = Ea - E;
   if dE < 0,               % Improved behavior
      ab = aux;
   elseif dE > 0,
      if rand < exp(-dE/T), % Probabilistic behavior
         % disp('** Annealing **');
         ab = aux;
      else,
         Eb = Ea;
      end;
   elseif rem(it,5) == 0,   % Identify Steady-State
      % disp('** Cooling **');
      T = alfa * T;
   end;
   if rem(it,20) == 0,
      disp(sprintf('It: %d	T: %2.3f	LR: %2.3f	E: %2.4f	Rb: %2.2f',it,T,lr,E,Rb));
      if lr > 5e-3,
         lr = 0.9*lr;
      end;
   end;
   if rem(it,50) == 0,
      figure(2); plotvet(ab); plot([0,0],[0,0],'ro'); title('Unit Directions'); drawnow;
   end;
   it = it + 1; Ea = E; vE = [vE;E]; vRb = [vRb;Rb];
end;
disp(sprintf('It: %d	T: %2.3f	LR: %2.3f	E: %2.4f	Rb: %2.2f',it,T,lr,E,Rb));
figure(2); plotvet(ab); plot([0,0],[0,0],'ro'); title('Unit Directions'); drawnow;

% --------------------- %
% INTERNAL SUBFUNCTIONS
% --------------------- %

% Function inserts noise into the vectors 
% & Keep them with unit norm
function [ab] = noise(ab,lr);
[N,s] = size(ab);
mat = rand(N,s) <= lr; % positions where to insert noise
aux = 10.*randn(N,s);
ab = ab + lr .* mat .* aux;
ab = unit_norm(ab);

% Function performs the unit vector normalization
function [ab] = unit_norm(ab);
[N,s] = size(ab);
norma = sqrt(sum(ab'.^2));
ab = ab ./ (ones(s,1)*norma)';

% Function determines the energy value
function [E,Rb] = ener(ab);
[N,s] = size(ab);
Ib = (1/N)*sum(ab);  % Average vector
Rb = sqrt(Ib*Ib');   % Norm of the average vector Ib
F = sum(sum(dist(ab,ab')));
Rb = 100 * (1 - Rb);
% E = (1/Rb) + F;
E = F;

% Function plots 2D vectors
function [] = plotvet(ab);
[N,s] = size(ab);
if s == 2,
   x = zeros(1,2); y = zeros(1,2);
   clf; xlabel('x'); ylabel('y'); % title('Unit directions');
   for i=1:N,
      x(1,2) = ab(i,1);	% coord. em x
      y(1,2) = ab(i,2);	% coord. em y
      hold on; plot(x,y);
      % plot(aux1,aux2,'co');
   end;
   axis([-1.1 1.1 -1.1 1.1]); 
elseif s == 3,
   x = zeros(1,2); y = zeros(1,2); z = zeros(1,2);
   clf;
   for i=1:N,
      x(1,2) = ab(i,1);	% coord. em x
      y(1,2) = ab(i,2);	% coord. em y
      z(1,2) = ab(i,3);	% coord. em z
      plot3(x,y,z); hold on;
   end;
   xlabel('x'); ylabel('y'); zlabel('z'); % title('Unit directions');
   axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]); 
end;
