function [rhoP1,rhoM1] = mexTest
cover = single(imread('E:\astego\Images\test\stegos\195.bmp'));
% t0 = tic;
[rhoP1,rhoM1] = CostHUGO(cover, single([1.1,1.5]));
% disp(toc(t0));