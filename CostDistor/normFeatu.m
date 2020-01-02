% function [rhoP1,rhoM1] = mexTest
%{
% E:\astego\Images\standard_images\covers\
root = 'E:\astego\Images\BOSS_ALL\';
name = '55.pgm';
cpath = [root,name];
SUNWDPath = ['E:\astego\StandExpers\SUNWD_04\',name];
HILLPpath = ['E:\astego\StandExpers\HILL_04\',name];
payLoad = single(0.4);

%% 隐写
cover = single(imread([root,name]));
stego = S_UNIWARD(uint8(cover), payLoad);
imwrite(uint8(stego), SUNWDPath, 'pgm');
stego = HILL(cover, payLoad);
imwrite(uint8(stego), HILLPpath, 'pgm');

%% 提取特征
Fc= SRMProces(SRM({cpath}), 0);
S_HILL= SRMProces(SRM({HILLPpath}),0);
S_SUNWD= SRMProces(SRM({SUNWDPath}),0);
DHILL= norm(S_HILL-Fc);
DSUNWD= norm(S_SUNWD-Fc);
%}
%% 失真度量
D = S_SUNWD04_SRM.F - C_Stand_SRM.F;
% D= getSRMQ1FromSRM(S_SUNWD_04_SRM.F) - getSRMQ1FromSRM(C_StandImg_SRM.F);
D_SUNWD = zeros(size(D,1),1);
for i=1:size(D,1)
  D_SUNWD(i)=norm(D(i,:));
end
clear i t0 root name payLoad cpath ans stego cover D;

%%
%{
stego = HUGO_like(uint8(cover), payLoad);
stego = S_UNIWARD(uint8(cover), payLoad);
stego = HILL(cover, payLoad);

[rhoP1,rhoM1] = CostHILL(cover);
[stego,rhoP1]=embedAlgCZL(cover, payLoad);
P = 1./rhoP1;


H = single([-1,2,-1; 2,-4,2; -1,2,-1]);
F = [-1, 2, -2, 2,-1;
      2,-6,  8,-6, 2;
     -2, 8,-12, 8,-2;
      2,-6,  8,-6, 2;
     -1, 2, -2, 2,-1];

% HF = filter2(F, I, 'same');
HF = imfilter(I,F);
figure;imshow(HF,[]);title('HUGO修改概率');

% [rhoP1,rhoM1] = CostHILL(I);
% P = 1./rhoP1;

D= getSRMQ1FromSRM(S_HILL_04_SRM.F) - getSRMQ1FromSRM(C_BOSS_SRM.F);
DHILL= zeros(size(D,1),1);
for i=1:size(D,1)
  DHILL(i)=norm(D(i,:));
end

[t, ~] = wiener2(I, [3, 3]);
HF = I - t;
%}