close all;

root = 'E:\astego\Images\tmp\stego\';
src1 = imread([root,'1.pgm']);
src2 = imread([root,'195.pgm']);
    
sHILL = HILL(srcPath, payLoad);
sUNWD = S_UNIWARD(uint8(src), payLoad);
sHUGO = HUGO_like(uint8(src), payLoad);

% resid
DHILL = sHILL - src;
DUNWD = single(sUNWD) - src;
DHUGO = single(sHUGO) - src;

figure('name','HILL'); imshow(DHILL,[]);
figure('name','UNWD'); imshow(DUNWD,[]);
figure('name','HUGO'); imshow(DHUGO,[]);