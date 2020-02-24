root = 'E:\Desktop\';
name = '1013.bmp'; cpath = [root,name];
cover= single(imread(cpath));

% correl(data);
% [rhoP1,~] = CostCZL(cover);
[stego,rhoP1]=CostCZLT11_mex(cover,9);

% P = 1./rhoP1;
% figure;histogram(P);
% figure; imshow(P,[]);