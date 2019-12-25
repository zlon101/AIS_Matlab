% 混合多种隐写算法
numSamp = size(C_BOSS_SRM.F,1);
b1 = 3000; b2 = 6000;
randSamps = randperm(numSamp);
S.names = cell(numSamp,1);
S.F= zeros(size(C_BOSS_SRM.F),'single');

inds = randSamps(1:b1);
S.names(1:b1) = S_HUGO_04_SRM.names(inds);
S.F(1:b1,:) = S_HUGO_04_SRM.F(inds, :);

inds = randSamps(b1+1:b2);
S.names(b1+1:b2) = S_SUNWD_04_SRM.names(inds);
S.F(b1+1:b2,:) = S_SUNWD_04_SRM.F(inds, :);

inds = randSamps(b2+1:end);
S.names(b2+1:end) = S_SUNWD_04_SRM.names(inds);
S.F(b2+1:end,:) = S_SUNWD_04_SRM.F(inds, :);

isequal(C_BOSS_SRM.names,sort(S.names))
[PE,PFA,PMD,trained_ensemble] = tutorial(C_BOSS_SRM, S);