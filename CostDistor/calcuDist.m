function [distortion,residual] =  calcuDist(coverImg,stegoImg)
% ÒþÐ´Ê§Õæ
% 
%%
% [rhoP1,rhoM1] = CostHUGO(coverImg);
% [rhoP1,rhoM1] = CostHILL(coverImg);
[rhoP1,rhoM1] = CostUNWD(coverImg);

residual = stegoImg - coverImg;
rhoM1 = rhoM1(residual==-1);
rhoP1 = rhoP1(residual==1);
distortion = sum(rhoM1) + sum(rhoP1);
end