function [bestFits,bestAbs,meanFits] = mexTest
payLoad = single(0.4);
Root = 'E:\astego\CSA\';
name = 'xx.pgm';
srcPath= [Root,'stegos\', name];
sharpedPath = [Root,'sharpeds\',name];
sharpedStegoPath = [Root,'sharpedStegos\',name];
embedParas = struct('srcPath',srcPath,'sharpedPath',sharpedPath,...
    'sharpedStegoPath',sharpedStegoPath,'payLoad',payLoad);
Memory.K = {};  Memory.V = zeros(20,1,'single');
Memory.last=uint8(1);
  
% [bestFits,bestAbs,meanFits] = CSA('E:\astego\Images\BOSS_ALL\1.pgm',...
%   single(0.4));
Abs = zeros(2,1,'single');

[fits,Mem]=calcuFit(Abs,embedParas,Memory);