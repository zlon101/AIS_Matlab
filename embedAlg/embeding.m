root='E:\astego\Images\BOSS_ALL\';
% root = 'E:\Desktop\';
name = '1013.pgm'; cpath = [root,name];
payload=single(0.4);
cover= single(imread(cpath));

% [rhoP1,rhoM1,optP1,optM1] = CostCZL(cover);
stego=embedAlgCZL(cpath,payload);
% diff=stego-cover;

% [rhoP1,rhoM1,optP1,optM1] = CostCZL(cover);
% P = 1./rhoP1;
% figure;histogram(P);
% figure; imshow(P,[]);