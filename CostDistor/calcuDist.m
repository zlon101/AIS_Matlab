function [distortion,residual] =  calcuDist(coverPath,stegoPath)
% ¼ÆËã stego-cover µÄÊ§Õæ
% 
%%
params.gamma = 1;
params.sigma = 1;
% [rhoP1,rhoM1] = CostHUGO(coverPath, params);
% [rhoP1,rhoM1] = CostUNIWD(coverPath);
[rhoP1,rhoM1] = CostHILL(coverPath);

residual = int8(imread(stegoPath)) - int8(imread(coverPath));
distM1 = rhoM1(residual==-1);
distP1 = rhoP1(residual==1);
distortion = sum(distM1) + sum(distP1);
end