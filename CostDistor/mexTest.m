function [rhoP1,rhoM1] = mexTest
cover = single(imread('E:\astego\Images\test\195.pgm'));
[rhoP1,rhoM1] = CostHUGO(cover, single([1.1,1.5]));