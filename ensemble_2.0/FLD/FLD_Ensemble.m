function learner = FLD_Ensemble(C, S, learner)
if(isstruct(C) && isstruct(S))
    C = C.F;    S = S.F;
end
% OOB.SUB = floor(size(C,1)*rand(size(C,1), 1)) + 1;
if( ~exist('learner','var') )
    learner.subspace = randperm(size(C,2));
end
Xm = C(:, learner.subspace);
Xp = S(:, learner.subspace);
% Xm = C;  Xp = S;
clear C S;

% remove constants
remove = false(1,size(Xm,2));
% 去除重复项
adepts = unique([find(Xm(1,:)==Xm(2,:)) find(Xp(1,:)==Xp(2,:))]);
for ad_id = adepts
    U1=unique(Xm(:,ad_id));
    if numel(U1)==1
        U2=unique(Xp(:,ad_id));
        if numel(U2)==1, if U1==U2, remove(ad_id) = true; end; end
    end
end
clear adepts U1 U2;

muC  = sum(Xm,1); muC = double(muC)/size(Xm,1);
muS  = sum(Xp,1); muS = double(muS)/size(Xp,1);
mu = (muS-muC)';

% calculate sigC
xc = bsxfun(@minus,Xm,muC);
sigC = xc'*xc;
% sigC = double(sigC)/size(Xm,1);  为了减低内存修改
sigC = single(sigC)/size(Xm,1);

% calculate sigS
xc = bsxfun(@minus,Xp,muS);
sigS = xc'*xc;
% sigS = double(sigS)/size(Xp,1);  为了减低内存修改
sigS = single(sigS)/size(Xp,1);
save('Xm', 'Xm');  save('Xp','Xp');
clear muC muS xc Xm Xp;

lengthOfSigC = length(sigC);
[rowsSigC, ~] = size(sigC);
sigCS = sigC + sigS;  % sigCS = sigC 34671*34671
clear sigS sigC;

% regularization
sigCS = sigCS + 1e-10*eye(rowsSigC);

% check for NaN values (may occur when the feature value is constant over images)
nan_values = sum(isnan(sigCS))>0;
nan_values = nan_values | remove;

sigCS = sigCS(~nan_values,~nan_values);
mu = mu(~nan_values);
lastwarn('');
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

% 投影矩阵 w=inv(sigCS * mu)
learner.w = sigCS\mu;
% regularization (if necessary)
[txt,warnid] = lastwarn(); %#ok<ASGLU>
while( strcmp(warnid,'MATLAB:singularMatrix') || (strcmp(warnid,'MATLAB:nearlySingularMatrix') ...
        && ~settings.ignore_nearly_singular_matrix_warning) )
    fprintf('\n\n-----ensemble_training.m 进入while循环------\n\n');
    % pause;
    lastwarn('');
    if ~exist('counter','var'), counter=1; else counter = counter*5; end
    sigCS = sigCS + counter*eps*eye(size(sigCS,1));
    learner.w = sigCS\mu;
    [txt,warnid] = lastwarn(); %#ok<ASGLU>
end    
warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:singularMatrix');
if( length(sigCS)~=lengthOfSigC )
    fprintf('\n\n-----ensemble_training.m length(sigCS)~=length(sigC)------\n\n');
    % pause;
    % resolve previously found NaN values, set the corresponding elements of w equal to zero
    w_new = zeros(lengthOfSigC,1);
    w_new(~nan_values) = learner.w;
    learner.w = w_new;
end
clear sigC mu;

% find threshold to minimize FA+MD
load('Xm');  load('Xp');
learner = findThreshold(Xm, Xp, learner);
learner.b = round(learner.b,2);

%{
figure;
P1 = Xm * learner.w;
P2 = Xp * learner.w;
vmean = mean(P1);
 border = learner.b;
% P1 = P1/vmean;  P2 = P2/vmean;  border = learner.b/vmean;
vmin = min(min(P1),min(P2));
vmax = max(max(P1),max(P2));
axis([0,length(P1), vmin, vmax]);
plot(P1, 1:length(P1), '*k');hold on;
plot(P2, 1:length(P2), 'xb');hold on;
hp = plot([border, border],[0,length(P1)], '-r');
% legend(hp,{['b: ', num2str(learner.b)]}, 'FontSize',12);
legend('cover','stego', ['b: ', num2str(learner.b)]);
%}
end

function [learner] = findThreshold(Xm,Xp,learner)
% find threshold through minimizing (MD+FA)/2, where MD stands for the
% missed detection rate and FA for the false alarms rate
P1 = Xm*learner.w;
P2 = Xp*learner.w;
L = [-ones(size(Xm,1),1);ones(size(Xp,1),1)];
[P,IX] = sort([P1;P2]);
L = L(IX);
Lm = (L==-1);
sgn = 1;

MD = 0;
FA = sum(Lm);
MD2=FA;
FA2=MD;
Emin = (FA+MD);
Eact = zeros(size(L-1));
Eact2 = Eact;
for idTr=1:length(P)-1
    if L(idTr)==-1
        FA=FA-1;
        MD2=MD2+1;
    else
        FA2=FA2-1;
        MD=MD+1;
    end
    Eact(idTr) = FA+MD;
    Eact2(idTr) = FA2+MD2;
    if Eact(idTr)<Emin
        Emin = Eact(idTr);
        iopt = idTr;
        sgn=1;
    end
    if Eact2(idTr)<Emin
        Emin = Eact2(idTr);
        iopt = idTr;
        sgn=-1;
    end
end

% 类边界border
learner.b = sgn*0.5*(P(iopt)+P(iopt+1));
if (sgn==-1)
    learner.w = -learner.w; 
end
end

