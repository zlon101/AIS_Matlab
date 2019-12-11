function [rhoP1,rhoM1] = mexTest
cover = single(imread('E:\astego\Images\test\stegos\195.bmp'));
[rhoP1,rhoM1] = CostHUGO(cover, single([1.1,1.5]));






%{
% create mirror padded cover image
% padSize = 3;
coverPadded = padarray(cover, [3,3], 'symmetric');
% create residuals
C_Rez_H = coverPadded(:, 1:end-1) - coverPadded(:, 2:end);
C_Rez_V = coverPadded(1:end-1, :) - coverPadded(2:end, :);
C_Rez_Diag = coverPadded(1:end-1, 1:end-1) - coverPadded(2:end, 2:end);
C_Rez_MDiag = coverPadded(1:end-1, 2:end) - coverPadded(2:end, 1:end-1);
stegoPadded = coverPadded;
% create residuals
S_Rez_H = stegoPadded(:, 1:end-1) - stegoPadded(:, 2:end);
S_Rez_V = stegoPadded(1:end-1, :) - stegoPadded(2:end, :);
S_Rez_Diag = stegoPadded(1:end-1, 1:end-1) - stegoPadded(2:end, 2:end);
S_Rez_MDiag = stegoPadded(1:end-1, 2:end) - stegoPadded(2:end, 1:end-1);


[rhoP1,rhoM1]=fn(cover, C_Rez_H,C_Rez_V,C_Rez_Diag, C_Rez_MDiag,...
  S_Rez_H, S_Rez_V,S_Rez_Diag,S_Rez_MDiag,single([1.1,1.5]));
%}