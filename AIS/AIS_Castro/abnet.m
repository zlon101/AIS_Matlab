function [w,win,cwin,D] = abnet(ag,eps,comp,alfa,beta,pc,pm),

%
% Ph.D. Thesis
% Leandro Nunes de Castro
% June, 1999
% Pattern Recognition in the Immune System using a Growing SOM
% Bipolar Splitting/Pruning Self-Organizing Feature Map (GSOM)
% with Evolutionary Phase
% Main features: bipolar weights, Hamming Distance, Winner takes all
%                PHASE I: Growing followed by Pruning
%                PHASE II: Supervised Evolution
%
% function [w,win,cwin,D] = hybrid(ag,eps,comp,alfa,beta,pc,pm),
% w        -> weight matrix (Ab population)
% win      -> winner for each Ag (v)
% cwin     -> amount of winning of each individual (tau)
% D        -> hamming distance of each Ag with relation to its mapped class
% ag       -> antigen population to be recognized (n2xs2)
% eps      -> ball of stimulation
% comp     -> comparison: 1 for comparing complementary chains
%                         0 for comparing identical chains (Hamm. dist.)
% alfa     -> amount of bits to be changed
% beta     -> number of iterations for reducing the learning rate
%
% Auxiliar functions: COVER, UPDATE, SPLIT, PRUNE, MATCH, CADEIA, TESTGSOM
% The columns of w must be similar to each Ag
%

if nargin == 2,
   [n2,s2] = size(ag);
   comp = 0;
   alfa = 3;
   beta = 3;
   pc = 0.6;
   pm = 0.1;
end;

% Network parameters
ep = 0; alfa0 = alfa; TD = 1;
[np,ni] = size(ag); no = 1; vep = [0];
[C,maxno] = cover(ni,eps); vno = [1:1:no];
disp(sprintf('Coverage of each Ab: %d',C));
disp(sprintf('Initial number of classes: %d',no));
disp(sprintf('Possible number of classes: %d',maxno));
if maxno > np,
   maxno = np; disp(sprintf('Maximum number of classes (N): %d',np));
end;
% disp(sprintf('Affinity threshold: %d',eps));
disp(sprintf('Press any key to continue...'));
pause;
[w] = cadeia(ni,no,0,0,1);
max_ep = (beta + 1) * maxno;

% Network Definition
while (ep < max_ep & TD > 0)% & no < maxno),
   cwin = zeros(1,no); k = 0;
   vet = randperm(np);			% Assincronous
   while k < np,
      k = k+1; i = vet(k); D = [];
      [D,mXOR] = match(w',ag(i,:),comp);
      [v(k),ind] = min(D);
      cwin(ind) = cwin(ind) + 1;
      win(i) = ind;
      w = update(w,ind,alfa,mXOR(ind,:)');
   end;
   TD = sum(v);
   ep = ep + 1;
   % Growing Phase
   if (rem(ep,beta)==0),
      [w,no,alfa] = split(cwin,win,w,ag,eps,alfa,alfa0);
      vno = [vno no]; vep = [vep ep];
   end;
   % Pruning Phase
   [aux,indmin] = min(cwin);
   if aux == 0,
      [w,no,alfa] = prune(w,indmin,alfa0);
      vno = [vno no];
   end;
   % Learning rate decreasing
   if (ep > 0.05*max_ep & rem(ep,0.05*max_ep)==0),
      if alfa > 1,
         alfa = alfa - 1;
      end;
   end;
   disp(sprintf('IT: %4.0d  no: %d  LR: %d	TD: %d',ep,no,alfa,TD));
end;
[v,win,cwin,perc] = testgsom(w,ag,eps);
disp(sprintf('Percentage of misclassified Ag: %3.2f%%',perc));
disp('Minimal Antigenic Affinity (HD)'); disp(v);
disp('Concentration Level: '); disp(cwin);
disp(sprintf('Final Architecture: [%d,%d].',ni,no));
figure(1); plot(vep,vno); hold on; plot(vep,vno,'or'); axis([0 ep+1 0 no+1]);
title('Growing Evolution');xlabel('Iteration'); hold off;

% --------------------------- %
%    INTERNAL SUBFUNCTIONS    %
% --------------------------- %

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


% Function SPLIT
function [w,no,alfa] = split(cwin,win,w,ag,eps,alfa,alfa0)
[ni,no] = size(w);
[ind] = find(cwin > 1);  % which outputs map more than one Ag
if ~isempty(ind),
   [val,out] = max(cwin);
%   out = ind(1);
   v = find(win==out);
   Mag = ag(v,:);        % matrix of ag mapped in the same output
   D = match(Mag,w(:,out)',0);
   [aux,new] = max(D);
   if aux > eps,
      disp('** Growing **');
      if out == 1,
         w = [Mag(new,:)',w];
      elseif out == no,
         w = [w,Mag(new,:)'];
      else,
         w = [w(:,1:out),Mag(new,:)',w(:,out+1:end)];
      end;
      no = no + 1;
      alfa = alfa0;
   end;
end;
% End Function SPLIT


% Function TESTGSOM
function [v,win,cwin,k] = testgsom(w,ag,eps),
% disp('** Running the trained network **');
[np,ni] = size(ag); k = 0;
cwin = zeros(1,size(w,2));
for i=1:np,
   [D] = match(w',ag(i,:),0);
   [v(i),ind] = min(D);
   win(i) = ind;
   cwin(ind) = cwin(ind) + 1;
end;
k = 100 * (sum(v > eps) / np);
% End Function TESTGSOM


% Function PRUNE
function [w,no,alfa] = prune(w,ind,alfa0),
[ni,no] = size(w);
disp('** Pruning **');
if ind == 1,
   w = w(:,2:no);
elseif ind == no,
   w = w(:,1:no-1);
else,
   w = [w(:,1:ind-1) w(:,ind+1:no)];
end;
no = no - 1;
alfa = alfa0;
% End Function PRUNE

% Function COVER
function [C,no,eps] = cover(len,eps),
fat = fatorial(len);
C = 0;
while eps > len,
   disp(sprintf('Ball of stimulation bigger than chain length %d',len));
   eps = input('Enter a new ball of stimulation: ');
end;

for i=0:eps,
   C = C + (fat/(fatorial(i) * fatorial(len-i)));
end;
no = ceil((2^len)/C);
% End Function COVER


% Function FATORIAL
function fat = fatorial(m);
if m == 0,
   fat = 1;
elseif m < 0,
   disp('Negative value');
else,
   fat = prod(1:1:m);
end;
% End Function FATORIAL


% Function UPDATE
function [w] = update(w,ind,alfa,vXOR),
[ni,no] = size(w);
for j = 1:alfa,
   [val,pto] = max(vXOR);
   if val == 0,
      break;     % exit loop if vectors are equal
   end;
   w(pto,ind) = -1 * w(pto,ind);
   vXOR(pto) = 0;
end;
% End Function UPDATE


% Function MATCH
function [ms,mXOR] = match(ab,ag,comp)
if nargin == 2,
   comp = 0;                     % Hamming distance
end;
msc = [];   % ms complement
% Converting bipolar (-1,+1) strings to binary ones (0,+1)
ab = hardlim(ab);
ag = hardlim(ag);
% Using the XOR operator for calculating the match score
[n1,s1] = size(ab);
ag = ones(n1,1) * ag;					% Multiply the Antigen
mXOR = xor(ab,ag);
ms = sum(mXOR');
msc = 1 - ms;
if comp == 1,
   ms = msc;
end;
% End Function MATCH