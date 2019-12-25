root = 'E:\astego\Images\standard_images\covers\';
name = '1.pgm';
cpath = [root,name];
SUNWDPath = ['E:\astego\标准图像集实验\SUNWD_04\',name];
HILLPpath = ['E:\astego\标准图像集实验\HILL_04\',name];
payLoad = single(0.4);

% cover = single(imread([root,name]));

%% 隐写
% stego = S_UNIWARD(uint8(cover), payLoad);
% imwrite(uint8(stego), SUNWDPath, 'pgm');
% stego = HILL(cover, payLoad);
% imwrite(uint8(stego), HILLPpath, 'pgm');

%% 提取特征
t0=tic;
Fc = SRM({cpath}); fprintf('SRM耗时：');disp(toc(t0));
S_HILL=SRM({HILLPpath});
S_SUNWD=SRM({SUNWDPath});

clear t0 root name payLoad;

%% 失真度量
% D = S_HILL_04_SRM.F - C_StandImg_SRM.F;
% % D= getSRMQ1FromSRM(S_SUNWD_04_SRM.F) - getSRMQ1FromSRM(C_StandImg_SRM.F);
% D_HILL = zeros(size(D,1),1);
% for i=1:size(D,1)
%   D_HILL(i)=norm(D(i,:));
% end

%%
%{
stego = HUGO_like(uint8(cover), payLoad);
stego = S_UNIWARD(uint8(cover), payLoad);
stego = HILL(cover, payLoad);

[rhoP1,rhoM1] = CostHILL(cover);
[stego,rhoP1]=embedAlgCZL(cover, payLoad);
P = 1./rhoP1;
figure;histogram(P);
%}