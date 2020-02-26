%% 1000幅载体图像作为图像库
n= 1000;
C= C_BOSS_SRM; S= Feat;
C2.names= C.names(1:n); C2.F= C.F(1:n,:);
S2.names= S.names(1:n); S2.F= S.F(1:n,:);

PE = tutorial(C2, S2);
clearvars -except PE C_BOSS_SRM Feat;

%{
root='D:\Matlab2017b\bin\featureDatas\';
% load([root,'covers\C_BOSS_SRM.mat']);
% load([root,'S_MiPOD_04_SRM.mat']); 
% load([root,'stegos\S_WOW_04_SRM.mat']);
C = C_BOSS_SRM;
S1 = S_BestArgsHUGO_04_SRM;
S2 = S_HUGO_04_SRM;
%% 混合多种隐写算法
% S_MiPOD_04_SRM S_WOW_04_SRM
numSamp = size(C.F,1);
b1 = 5000; b2 = 6000;
randSamps = randperm(numSamp);
S.names = cell(numSamp,1);
S.F= zeros(size(C.F),'single');

inds = randSamps(1:b1);
S.names(1:b1) = S1.names(inds);
S.F(1:b1,:) = S1.F(inds, :);

inds = randSamps(b1+1:end);
S.names(b1+1:end) = S2.names(inds);
S.F(b1+1:end,:) = S2.F(inds, :);

% inds = randSamps(b2+1:end);
% S.names(b2+1:end) = S_SUNWD_04_SRM.names(inds);
% S.F(b2+1:end,:) = S_SUNWD_04_SRM.F(inds, :);

clearvars -except C S;
% fprintf('isequal: '); disp(isequal(C.names,sort(S.names)));
[PE,PFA,PMD] = tutorial(C, S);
%}