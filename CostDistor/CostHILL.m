function [rhoP1,rhoM1] = CostHILL(coverPath)
%COSTUNIWD 此处显示有关此函数的摘要
%   此处显示详细说明
%%
src = single(imread(coverPath));
[M,N] = size(src);
wetCost = 10^8;
%% 高通H-KB滤波器
H = single([-1,2,-1; 2,-4,2; -1,2,-1]);
Y = filter2(H, src, 'same');
%% 低通L1
L1 = single(ones(3));
Z = filter2(L1, abs(Y), 'same');
%% 低通L2
L2 = single(ones(15));
Cost = filter2(L2, Z.^-1, 'same');
%% adjust embedding costs
Cost(Cost > wetCost) = wetCost; % threshold on the costs
Cost(isnan(Cost)) = wetCost; % if all xi{} are zero threshold the cost
rhoP1 = Cost;  % +1 的代价
rhoM1 = Cost;
rhoP1(src==255) = wetCost; % do not embed +1 if the pixel has max value
rhoM1(src==0) = wetCost; % do not embed -1 if the pixel has min value
end