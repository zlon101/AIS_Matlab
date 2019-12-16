function [D,R] =  calcuDist(coverImg,stegoImg)
% ÒþÐ´Ê§Õæ
% 
%%
t0= tic;
Fc= SRMProces(SRMQ1({coverImg}),0);
disp(toc(t0));
Fs= SRMProces(SRMQ1({stegoImg}),0);
D = norm(Fs-Fc);
R = 0;

% [rhoP1,rhoM1] = CostHUGO(coverImg);
% [rhoP1,rhoM1] = CostUNWD(coverImg);
% [rhoP1,rhoM1] = CostHILL(coverImg);
% R = stegoImg - coverImg;
% rhoM1 = rhoM1(R==-1);
% rhoP1 = rhoP1(R==1);
% D = sum(rhoM1) + sum(rhoP1);
end