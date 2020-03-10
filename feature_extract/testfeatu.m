% -------------------------------------------------------------------------
% I = imread('E:\astego\Images\BOSS_ALL\1.pgm');
% S = ones(size(I));
% t0= datetime('now');
% F = maxSRM(I, S);
% disp(datetime('now')-t0);
% -------------------------------------------------------------------------
root= 'E:\astego\Images\BOSS_ALL\';
outPath='E:\astego\CZL\F.mat';
maxSRMEXE(root,outPath,'1','3')