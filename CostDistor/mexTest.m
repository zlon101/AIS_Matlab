% function [rhoP1,rhoM1] = mexTest
% close all;
root = 'E:\astego\Images\Experis\';
name = '195.bmp';
cPath = [root,'1covers\',name];
sPath = [root,'HILL04\',name];
cover = single(imread(cPath));
[rhoP1,rhoM1] = CostHUGO(cover);

P = 1./rhoM1;
figure('name','HUGO');histogram(P);
figure('name','HUGO');imshow(P,[]);
clear root name cPath sPath cover;
% figure('name','µÈÈ«_-1¸ÅÂÊ');imshow(P,[]);