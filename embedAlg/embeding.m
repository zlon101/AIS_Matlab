root='E:\astego\Images\BOSS_ALL\';
name = '1013.pgm'; cpath = [root,name];
payload=single(0.4);
cover= single(imread(cpath));

% stego=embedAlgCZL(cpath,payload);
% [rhoP1,rhoM1] = CostCZL(cover);

%--------------------------------------
sigma= 1;
Fsize= 2*ceil(2*sigma)+1;

Filter1 = fspecial('gaussian',Fsize,sigma);
Filter2 = fspecial('gaussian',3,sigma);
%--------------------------------------