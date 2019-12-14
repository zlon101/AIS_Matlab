function [bestFits,bestAbs,meanFits] = mexTest
payLoad = single(0.4);
Root = 'E:\astego\Images\test\';
name = '195.pgm';
stego = HUGO_like(imread([Root,name]),payLoad);
imwrite(uint8(stego),[Root,'195_HUGOLike.pgm'],'pgm');

srcPath= [Root,'stegos\', name];
sharpedPath = [Root,'sharpeds\',name];
sharpedStegoPath = [Root,'sharpedStegos\',name];
embedParas = struct('srcPath',srcPath,'sharpedPath',sharpedPath,...
    'sharpedStegoPath',sharpedStegoPath,'payLoad',payLoad);
Memory.K = {};  Memory.V = zeros(20,1,'single');
Memory.last=uint8(1);
Abs = zeros(2,1,'single');
[fits,Mem]=calcuFit(Abs,embedParas,Memory);