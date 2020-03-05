root='E:\astego\Images\BOSS_ALL\';
name = '1013.pgm'; cpath = [root,name];
payload=single(0.4);
cover= single(imread(cpath));

% [rhoP1,rhoM1] = CostHUGO_like(cover);
[rhoP1,rhoM1] = CostCZL(cover);
% stego=embedAlgCZL(cpath,payload);

%---------------------------------------------------------------------------
% rhoP1= rhoP1(154:291,58:228);
P2=1./rhoP1;
% figure('name','HUGO');histogram(P1);
% figure('name','HUGO');imshow(P1,'DisplayRange',[0,1.5],'Border','tight');
% figure('name','CZL');histogram(P2);
figure('name','CZL');imshow(P2,'DisplayRange',[1,50],'Border','tight');
%---------------------------------------------------------------------------
% sigma= 13;
% Fsize= 2*ceil(2*sigma)+1;
% Fsize=13;
% Filter1 = fspecial('gaussian',Fsize,sigma);
%---------------------------------------------------------------------------