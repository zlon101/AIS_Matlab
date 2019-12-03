function [distortion,residual] =  calcuDist(coverImg,stegoImg)
% ÒþÐ´Ê§Õæ
% 
%%
% params.gamma = 1;  params.sigma = 1;
% [rhoP1,rhoM1] = CostHUGO(coverPath, params);
% [rhoP1,rhoM1] = CostUNIWD(coverPath);
[rhoP1,rhoM1] = CostHILL(coverImg);

% residual = single(stegoImg) - single(coverImg);
rhoM1 = rhoM1(single(stegoImg)-single(coverImg)==-1);
rhoP1 = rhoP1(single(stegoImg)-single(coverImg)==1);
distortion = sum(rhoM1) + sum(rhoP1);
residual = 0;
end