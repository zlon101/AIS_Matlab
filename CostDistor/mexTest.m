function [distortion,residual] = mexTest
% Matlab Coder ���� MEX�ļ�ʱʹ�õĲ����ļ�
root = 'E:\astego\Images\BOSS_ALL\';
src = single(imread([root,'1.pgm']));
stego = single(imread([root,'2.pgm']));
[distortion,residual] =  calcuDist(src,stego);