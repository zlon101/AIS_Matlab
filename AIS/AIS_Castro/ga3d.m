function [x,y,fx,vx,vmfit,P] = ga3d(v,ger,pc,pm);
%
% Ph.D. Thesis
% Leandro Nunes de Castro
% November, 1999.
% ENHANCED GENETIC ALGORITHM - Bi-classist Selection
%
% Secondary functions: DECODE & IMPRIME(internal)
%
% function [x,y,fx,vx,vmfit,P] = ga3d(v,ger,N,pc,pm);
% x,y,fx 	-> f(x,y) is the maximal value of the function
% vfx		-> vector of best fitness of population through generations
% vmfit	-> vector of mean fitness of population through generations
% v    	-> population of size N x 2.L (L = 22 [-1,1])
% ger		-> number of generations
% N			-> population size
% pc		-> cross-over probability
% pm		-> mutation probability
%
% Definitions
% f			-> evaluation function
% v			-> population
% fit		-> fitness vector
% x,y		-> decodified coordinates
% Maximum found: [x,y,f] = [1.76,1.27,4.0172];
%

% Default Parameters
if nargin == 1,
   [N,L] = size(v); ger = 25; pc = 0.5; pm = 0.01;
end;

disp(sprintf('Number of generations: %d',ger));
disp(sprintf('Population size: %d',N));
disp(sprintf('Crossover probability: %.3f',pc));
disp(sprintf('Mutation probability: %.3f',pm));
% f = '-1 * x .* sin(2 * pi .* x) + y.* sin(2 * pi .* y) + 1'; 
f = '1 * x .* sin(4 * pi .* x) - 1 * y.* sin(4 * pi .* y + pi) + 1'; 
% Plotting parameters
% [x,y] = meshgrid(-1:0.05:2,-1:0.05:2); vxp = x; vyp = y;
[x,y] = meshgrid(-1:0.05:1,-1:0.05:1); vxp = x; vyp = y;
vzp = eval(f); PRINT = 1;

% General parameters & Initial operations
sol = 1; vmfit = []; it = 1; vx = []; C = [];
x = decode(v(:,1:22)); y = decode(v(:,23:end)); fit = eval(f);
imprime(PRINT,vxp,vyp,vzp,x,y,fit,1,1); title('Initial Population');
disp('Press any key to continue...'); pause;
pb = 0.5; pw = 0.1; pr = 1 - (pb + pw);
nb = round(pb * N); nw = round(pw * N); nr = round(pr * N);
if (nb + nw + nr) ~= N,
   dif = N - (nb + nw + nr);
   nr = nr + dif;
end;

% Generations
t0 = clock; f0 = flops;
while it <= ger & sol <= 2.26,
	% 	Reproduction (Bi-classist Selection)
   [rw,ind] = sort(fit); ind = fliplr(ind);
   vtemp = [v(ind(1:nb),:); v(ind(end-nw+1:end),:); v(2:nr+1,:)];
   
   %  Crossover
   C(:,1) = rand(N,1) <= pc; C(:,2) = round(19.*rand(N,1))+1;
   I = find(C(:,1) == 1); IP = [I,C(I,2)];
   for i = 1:size(IP,1),
      v(IP(i,1),:) = [vtemp(IP(i,1),1:IP(i,2)) vtemp(1,IP(i,2)+1:end)];
   end;
         
   % 	Mutation
   M = rand(N,L) <= pm; M(1,:) = zeros(1,L);
   v = v - 2 .* (v.*M) + M;
   
   % Results
   % x = decode(v); fit = eval(f);
   x = decode(v(:,1:22)); y = decode(v(:,23:end)); fit = eval(f);
   imprime(PRINT,vxp,vyp,vzp,x,y,fit,it,5);
   [sol,indb] = max(fit);
	v(1,:) = v(indb,:); media = mean(fit);
	vx = [vx sol]; vmfit = [vmfit media];
	if rem(it,1) == 0 | it == 10,
		disp(sprintf('Gen.: %d  x: %2.2f  y: %2.2f  Av: %2.2f  f(x,y): %2.3f',it,x(indb),y(indb),media,sol));
	end;
	it = it + 1;
end;
T = etime(clock,t0); F = flops - f0;
x = x(indb); y = y(indb); fx = sol; P = v;

% 	Results
% clf; plot(vx); title('f(x) evolution'); xlabel('Generations'); ylabel('f(x)');
% hold on; plot(vmfit); hold off;

% --------------------- %
% INTERNAL SUBFUNCTIONS
% --------------------- %

% Print Surface
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

% Decodify bitstrings
function x = decode(v);
% x		-> real value (precision: 6)
% v		-> binary string (length: 22)
v = fliplr(v); s = size(v);
aux = 0:1:21; aux = ones(s(1),1)*aux;
x1 = sum((v.*2.^aux)');
x = -1 + x1 .* (2 / 4194303);
