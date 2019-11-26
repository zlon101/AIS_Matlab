function [x,y,fx,vfx,vmfit,P,vpm] = imalg3d(P,gen,n,pm,per)
% Ph.D. Thesis
% Leandro Nunes de Castro
% November, 1999.
% Immune Algorithm - New evolutionary strategy inspired in the Immune System
% Operations: Hypermutation, Editing, Selection
% Each clone has size proportional to its affinity
% IMALG3D: MULTI-PEAK Solution (3D Function)
%
% function [x,y,fx,vfx,vmfit,P] = imalg3d(P,gen,n,pm,per);
% x,y,fx 	-> f(x,y) is the maximal value of the function
% vfx		-> vector of best fitness of population through generations
% vmfit	-> vector of mean fitness of population through generations
% P    	-> population of size N x 2.L
% gen  	-> generation number
% n    	-> number of clones
% pm   	-> hypermutation probability
% per  	-> percentile of the population to suffer random reshuffle
% 
% T    	-> temporary population
% Maximum found: [x,y,f] = [1.76,1.27,4.0172];
%
if nargin == 1
   % gen = 200; n = round(size(P,1)/2); pm = 0.0005; per = 0.0; fat = 10;
   gen = 25; n = size(P,1); pm = 0.01; per = 0.0; fat = .1;
end
while n <= 0
   n = input('n has to be at least one. Type a new value for n: ');
end

f = '1 * x .* sin(4 * pi .* x) - 1 * y.* sin(4 * pi .* y + pi) + 1';
[x,y] = meshgrid(-1:0.05:1,-1:0.05:1); vxp = x; vyp = y;
vzp = eval(f);
x = decode(P(:,1:22)); y = decode(P(:,23:end)); fit = eval(f);
imprime(1,vxp,vyp,vzp,x,y,fit,1,1); title('Initial Population');
% disp(sprintf('Number of generations: %d',gen));
% disp(sprintf('Population size: %d',n));
% disp(sprintf('Mutation probability: %.3f',pm));
% disp(sprintf('Number of clones per candidate: %d',fat*n));
% disp(sprintf('Press any key to continue...')); 
% pause;

% Hypermutation controlling parameters
pma = pm; itpm = gen; pmr = 0.8;
% General defintions
vpm = []; vfx = []; vmfit = []; valfx = 1; 
[N,L] = size(P); it = 0; PRINT = 1;

% Generations
while it <= gen && valfx <= 2.26
   x = decode(P(:,1:22)); y = decode(P(:,23:end)); T = []; cs = [];
   fit = eval(f); [a,ind] = sort(fit);
   valx = x(ind(end-n+1:end)); valy = y(ind(end-n+1:end));
   fx = a(end-n+1:end);	% n best individuals (maximization)
   imprime(PRINT,vxp,vyp,vzp,x,y,fit,it,5);
   
   % Reproduction
   [T,pcs] = reprod(n,fat,N,ind,P,T);
   
   % Hypermutation
   M = rand(size(T,1),L) <= pm;
   T = T - 2 .* (T.*M) + M;
   T(pcs,:) = P(fliplr(ind(end-n+1:end)),:);
   
   % New Re-Selection (Multi-peak solution)
   x = decode(T(:,1:22)); y = decode(T(:,23:end)); fit = eval(f); pcs = [0 pcs];
   for i=1:n
      [out(i),bcs(i)] = max(fit(pcs(i)+1:pcs(i+1)));		% Maximizationn problem
      bcs(i) = bcs(i) + pcs(i);
   end
   P(fliplr(ind(end-n+1:end)),:) = T(bcs,:);
   
   % Editing (Repertoire shift)
   nedit = round(per*N); it = it + 1;
   P(ind(1:nedit),:) = cadeia(nedit,L,0,0,0);
   pm = pmcont(pm,pma,pmr,it,itpm); valfx = max(fx);
   vpm = [vpm pm]; vfx = [vfx valfx]; vmfit = [vmfit mean(fit)];
   disp(sprintf('It.: %d  pm: %.4f  x: %2.2f  y: %2.2f  Av.: %2.2f  f(x,y): %2.3f',it,pm,valx(end),valy(end),vmfit(end),valfx));
% end while
end
%imprime(PRINT,vxp,vyp,vzp,x,y,fit,it,1);
x = valx(end); y = valy(end); fx = max(fx);
% x = P(ind(end),1:22); y = P(ind(end),23:44); fx = max(fx);

% --------------------- %
% INTERNAL SUBFUNCTIONS
% --------------------- %

% Print
function [] = imprime(PRINT,vx,vy,vz,x,y,fx,it,mit);
% x,fx				-> actual values
% vxplot, vplot	-> original (base) function
if PRINT == 1,
   if rem(it,mit) == 0,
      mesh(vx,vy,vz); hold on; axis([-1 1 -1 1 -1 2.5]);
      xlabel('x'); ylabel('y'); zlabel('f(x,y)');
      plot3(x,y,fx,'k*'); drawnow; hold off;
   end;
end;

% Reproduction
function [T,pcs] = reprod(n,fat,N,ind,P,T);
% n		-> number of clones
% fat	-> multiplying factor
% ind	-> best individuals
% T		-> temporary population
% pcs	-> final position of each clone
if n == 1,
   cs = N;
   T = ones(N,1) * P(ind(1),:);
else,
   for i=1:n,
      % cs(i) = round(fat*N/i);
      cs(i) = round(fat*N);
      pcs(i) = sum(cs);
      T = [T; ones(cs(i),1) * P(ind(end-i+1),:)];
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

% Decodify bitstrings
function x = decode(v);
% x		-> real value (precision: 6)
% v		-> binary string (length: 22)
v = fliplr(v); s = size(v);
aux = 0:1:21; aux = ones(s(1),1)*aux;
x1 = sum((v.*2.^aux)');
x = -1 + x1 .* (2 / 4194303);

% Function CADEIA
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
% End Function CADEIA

   
