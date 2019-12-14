root = 'E:\astego\Images\Experis\';
name = '195.pgm';
cpath = [root,name];
spath = [root,'\HILL\',name];
payLoad = single(0.2);

cover = single(imread(cpath));
% [rhoP1,rhoM1] = CostHUGO(cover);


alg = 'HUGOOPT';
t0=tic;
stego = HUGO(cover, payLoad);
%{
disp(toc(t0));
% stego = HUGO_like(uint8(cover), payLoad);
% stego = S_UNIWARD(uint8(cover), payLoad);
% stego = HILL(cover, payLoad);
figure('name',alg);imshow(single(stego)-cover,[]);

[rhoP1,rhoM1] = CostHUGO(cover);
P = 1./rhoM1;
figure('name',alg);imshow(P,[]);

clear root name cpath spath payLoad alg;


imwrite(uint8(sHUGO), [root,'\HUGO\',name], 'pgm');
imwrite(uint8(sUNWD), [root,'\UNWD\',name], 'pgm');
imwrite(uint8(sHILL), [root,'\HILL\',name], 'pgm');

[DHUGO,RHUGO] = calcuDist(src, single(sHUGO));
[DUNWD,RUNWD] = calcuDist(src, single(sUNWD));
[DHILL,RHILL] = calcuDist(src, single(sHILL));
%}