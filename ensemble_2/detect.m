root='D:\Matlab2017b\bin\featureDatas\';
% load([root,'covers\C_BOSS_SRM.mat']);
load([root,'S_MiPOD_04_SRM.mat']); 
load([root,'stegos\S_WOW_04_SRM.mat']);

% 混合多种隐写算法
% S_MiPOD_04_SRM S_WOW_04_SRM
numSamp = size(C_BOSS_SRM.F,1);
b1 = 5000; b2 = 6000;
randSamps = randperm(numSamp);
S.names = cell(numSamp,1);
S.F= zeros(size(C_BOSS_SRM.F),'single');

inds = randSamps(1:b1);
S.names(1:b1) = S_MiPOD_04_SRM.names(inds);
S.F(1:b1,:) = S_MiPOD_04_SRM.F(inds, :);

inds = randSamps(b1+1:end);
S.names(b1+1:end) = S_WOW_04_SRM.names(inds);
S.F(b1+1:end,:) = S_WOW_04_SRM.F(inds, :);

% inds = randSamps(b2+1:end);
% S.names(b2+1:end) = S_SUNWD_04_SRM.names(inds);
% S.F(b2+1:end,:) = S_SUNWD_04_SRM.F(inds, :);

clearvars -except C_BOSS_SRM S;
isequal(C_BOSS_SRM.names,sort(S.names))
[PE,PFA,PMD] = tutorial(C_BOSS_SRM, S);