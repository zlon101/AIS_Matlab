function [] = tutorial
%
% Function TUTORIAL Demonstration
% Runs a Demo for the following immune tools:
% 1) ABNET (Growing Boolean Competitive Network - journal.doc)
% 2) CLONALG (Clonal Selection Algorithm - gecco00.doc)
% 3) AINET (Artificial Immune Network - ainet.doc)
%
% Other Sources of Reference: RT DCA 01/99 (Artificial Immune Systems Part I - Basic Theory and Applications)
%
% Secondary Functions: 
% Internal Subfunctions: PRESENT, CHOOSE, DRAW_ABNET, PLOTVET, CADEIA
% Auxiliary Subfunctions: ABNET, GA3D, IMALG3D, AINET, ANALYSIS, DENDRO
%
% Copyrigth by Leandro Nunes de Castro
% April/May, 2000
%

% Function chooses the demo
while 1,
   p = present; choose(p);
   disp('Press any key to continue...'); pause;
end;


% ------------------------- %
%   INTERNAL SUBFUNCTIONS
% ------------------------- %

% Function presents the DEMO TOOLBOX
function p = present,
clc; disp('*-----------------------------------------------*');
disp(sprintf('ARTIFICIAL IMMUNE SYSTEMS - TUTORIAL PRESENTATION'));
disp(sprintf('	    Leandro Nunes de Castro'));
disp(sprintf('		April/May, 2000'));
disp('*-----------------------------------------------*');
disp(sprintf('Available Demos:'));
disp(sprintf('1) SAND'));
disp(sprintf('2) ABNET'));
disp(sprintf('3) CLONALG'));
disp(sprintf('4) AINET'));
p = input('Type the desired demo number (or CTRL^C to Interrupt): ');
% End function PRESENT


% Function allows choosing Demos
function choose(p);
switch p,
case 1,
   disp(sprintf('\n** PART I - SAND (Simulated ANnealing Approach to Increase Diversity) **'));
   disp(sprintf('Available Demo Tasks:'));
   disp(sprintf('1) REAL-VALUED VECTORS'));
   disp(sprintf('2) 3-D BINARY POINTS'));
   task = input('Type the desired demo number (1) or (2): ');
   switch task,
   case 1,
      N = input('Choose the number of 2-D vectors: ');
      ab = randn(N,2);
      [ab,vE,vRb,E,Rb] = sareal(ab);
      figure(3); plot(vRb); title('Energy Evolution (E)');
      axis([0 length(vRb) 0 105]);
   case 2,
      N = input('Choose the population size (even size) [2,8]: ');
      [ab,vE,E] = sacubic(N);
      figure(3); plot(vE); title('Energy Evolution (E)');
      axis([0 length(vE) 0 105]);
   otherwise,
      disp('Accepted values are (1) or (2)');
   end; % End switch task SAD
case 2,
   disp(sprintf('\n** PART II - ABNET (An AntiBody NETwork) **'));
   L = 10; aux1 = ones(1,L); aux2 = zeros(1,L);
   ag = [aux1;[0,aux1(2:end)];[0,0,aux1(3:end)];[aux1(1:end-1),0];[aux1(1:end-2),0,0];...
         aux2;[1,aux2(2:end)];[1,1,aux2(3:end)];[aux2(1:end-1),1];[aux2(1:end-2),1,1]];
   disp(sprintf('\nAntigens (columnwise) to be Recognized: [%d,%d]',size(ag,1),L)); 
   disp(ag'); ag = ag - 0.1; ag = hardlims(ag);
   eps = input('Choose the Affinity Threshold [0,L]: ');
   [w,win,cwin,D] = abnet(ag,eps);
   figure(2); clf; draw_abnet(w);
case 3,
   disp(sprintf('\n** PART III - CLONALG (The CLONal Selection ALGorithm) **'));
   % figure(1); clf;
   v = cadeia(100,44,0,0,0);
   disp(sprintf('Available Demo Tasks:'));
   disp(sprintf('1) GA (GLOBAL SEARCH)'));
   disp(sprintf('2) CLONALG (MULTI-MODAL OPTIMIZATION)'));
   disp(sprintf('3) CHARACTER (KNOWLEDGE ACQUISITION)'));
   ga_cga = input('Type the desired demo number (1) ,(2) or (3): ');
   switch ga_cga,
   case 1,
      disp(sprintf('** Standard Genetic Algorithm - GA **'));
      figure(1); clf;
      [x,y,fx,vx,vmfit,P] = ga3d(v);
      disp(sprintf('Maximum found [x,y,f(x,y)]: [%.2f,%.2f,%.2f]',x,y,fx));
      figure(2); plot(vx); title('f(x,y) x Mean'); xlabel('Generations'); ylabel('f(x)');
      hold on; plot(vmfit,'r'); hold off;
   case 2,
      disp(sprintf('** Clonal Selection Algorithm - CLONALG **'));
      figure(3); clf; %v = v(1:50,:);
      [x,y,fx,vfx,vmfit,P,vpm] = imalg3d(v);
      save data x y fx vfx;
      disp(sprintf('Maximum found [x,y,f(x,y)]: [%.2f,%.2f,%.2f]',x,y,fx));
      figure(4); clf; plot(vfx); title('f(x,y) x Mean'); xlabel('Generations'); ylabel('f(x)');
      hold on; plot(vmfit,'r'); hold off;
   case 3,
      disp(sprintf('** CLONALG - Knowledge Acquisition**'));
      load num8_12x10; X = v(1,:); P = cadeia(10,120,0,0,0);
      M = imalgchar(P,X);
   otherwise,
      display('Accepted values are (1), (2) or (3)');
   end; % End Switch ga_cga
case 4,
   disp(sprintf('\n** PART IV - AINET (Artificial Immune NETwork) **'));
   disp(sprintf('Available Demo Tasks:'));
   disp(sprintf('1) LINEARLY SEPARABLE CLASSES'));
   disp(sprintf('2) 2-SPIRALS'));
   disp(sprintf('3) CHAINLINK'));
   task = input('Type the desired demo number (1) ,(2) or (3): ');
   switch task,
   case 1,
      load Ag; Ag = norma(Ag);
      figure(1); plotvet(Ag,10); title('Training Patterns');
      ts = input('Choose the suppression threshold (suggested: 0.1): ');
      [M,D] = ainet(Ag,ts,4,10,10,0.2,1);
      [E,Nc,mCe,mC,T] = anal1(M,D);
   case 2,
      figure(1); Ag = spir; title('Training Patterns');
      ts = input('Choose the suppression threshold (suggested: 0.07): ');
      [M,D] = ainet(Ag,ts,4,10,20,0.1,1);
      [E,Nc,mC,T] = anal2(M,D);
      disp(sprintf('Number of clusters: %d',Nc));
   case 3,
      figure(1); Ag = donuts1(500); title('Training Patterns');
      ts = input('Choose the suppression threshold (suggested: 0.2): ');
      [M,D] = ainet(Ag,ts,4,10,5,0.1,1);
      [E,Nc,mC,T] = anal3(M,D);
      disp(sprintf('Number of clusters: %d',Nc));
  otherwise,
      disp('Accepted values are (1), (2) or (3)');
   end; % End Switch task
      
   otherwise,
   disp(sprintf('Accepted values are (1), (2) or (3)'));
   break;
end;
% End Function CHOOSE


% Function DRAW_ABNET
function draw_abnet(w);
[N,L] = size(w);
r = [0.2.*ones(1,N);1:1:N]';
c = [1.2.*ones(1,L);1:1:L]';
if N > L, y = N; else, y = L; end;
plot(r(:,1),r(:,2),'sk'); hold on;
plot(c(:,1),c(:,2),'or');
axis([0 1.4 0 y+1]);
for i = 1:N,
   for j = 1:L,
      if w(i,j) == 1,
         line([r(i,1),c(j,1)],[r(i,2),c(j,2)]);
      end;
   end;
end;
title('ABNET (0 weights ommitted)');
hold off;
% End Function DRAW_ABNET


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


% Function Plot Vectors
function [] = plotvet(X,nec)
% nec    -> number of elements per class
clf; hold on;
xlabel('x'); ylabel('y');
plot(X(1:nec,1),X(1:nec,2),'bs');
plot(X(nec+1:2*nec,1),X(nec+1:2*nec,2),'ro');
plot(X(2*nec+1:3*nec,1),X(2*nec+1:3*nec,2),'g+');
plot(X(3*nec+1:4*nec,1),X(3*nec+1:4*nec,2),'kx');
plot(X(4*nec+1:5*nec,1),X(4*nec+1:5*nec,2),'mv');
hold off;
% End Function PLOTVET


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