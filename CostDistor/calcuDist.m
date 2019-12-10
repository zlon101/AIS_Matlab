function [distortion,residual] =  calcuDist(coverImg,stegoImg)
% ÒþÐ´Ê§Õæ
% 
%%
% [rhoP1,rhoM1] = CostHUGO(coverImg);
% [rhoP1,rhoM1] = CostHILL(coverImg);
[rhoP1,rhoM1] = CostUNIWD(coverImg);

rhoM1 = rhoM1(stegoImg-coverImg==-1);
rhoP1 = rhoP1(stegoImg-coverImg==1);
distortion = sum(rhoM1) + sum(rhoP1);
residual = stegoImg - coverImg;
end