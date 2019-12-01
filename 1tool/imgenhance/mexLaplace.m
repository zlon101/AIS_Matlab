function [dstImg,HF] = mexLaplace()
imgData = single(imread('E:\astego\CSA\stegos\xx.pgm'));
[dstImg,HF] =  laplace(imgData, single(0.5));
end