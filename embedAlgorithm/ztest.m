src = 'E:\astego\Images\standard_test_images\bmp\';
dst = 'E:\astego\Images\standard_test_images\pgm\';


close all;
root = 'E:\astego\Images\Experis\';
name = '1004.pgm';
srcPath = [root, name];
payLoad = single(0.2);
srcData = single(imread(srcPath));

sHUGO = HUGO_like(uint8(srcData), payLoad); % 1
sUNWD = S_UNIWARD(uint8(srcData), payLoad); % 2
sHILL = HILL(srcPath, payLoad); % 3.pgm

imwrite(uint8(sHUGO), [root,'\stegos\1.pgm'], 'pgm');
imwrite(uint8(sUNWD), [root,'\stegos\2.pgm'], 'pgm');
imwrite(uint8(sHILL), [root,'\stegos\3.pgm'], 'pgm');

[D1,resid1] =  calcuDist(srcPath, [root,'\stegos\1.pgm']);
[D2,resid2] =  calcuDist(srcPath, [root,'\stegos\2.pgm']);
[D3,resid3] =  calcuDist(srcPath, [root,'\stegos\3.pgm']);


% figure('name','HILL'); imshow(DHILL,[]);
% figure('name','UNWD'); imshow(DUNWD,[]);
% figure('name','HUGO'); imshow(DHUGO,[]);