root = 'E:\astego\Images\Experis\';
name = '195.pgm';
cpath = [root,name];
spath = [root,'\HILL\',name];
payLoad = single(0.4);

cover = single(imread(cpath));

[rhoP1,rhoM1] = CostHILL(cover);
% [stego,rhoP1]=embedAlgCZL(cover, payLoad);
P = 1./rhoP1;
figure;histogram(P);
% figure;imshow(P);

%{
disp(toc(t0));
% stego = HUGO_like(uint8(cover), payLoad);
% stego = S_UNIWARD(uint8(cover), payLoad);
% stego = HILL(cover, payLoad);
%}