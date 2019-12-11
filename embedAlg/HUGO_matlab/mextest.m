function mextest()
% matlab coder mex function
cover = single(imread('E:\astego\Images\test\stegos\195.bmp'));
[rhoP1,rhoM1] = CostHUGO(cover, single([1.1,1.5]));

