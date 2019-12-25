% This example demonstrates how to use the MiPOD embedding function
clc
clear all
close all

% Read the input cover image
Cover = double(imread ('1.pgm'));

% Set the payload to 0.4 bpp
Payload = 0.4;

% MiPOD embedding
tStart = tic;
[Stego, pChange, ChangeRate] = MiPOD( Cover, Payload );
tEnd = toc(tStart);
fprintf('MiPOD embedding is done in: %f (sec)\n',tEnd);

%%
close all

figure;
imshow (Cover,[]);
title ('Cover image');

figure;
imshow(1-pChange/0.3333);
title('MiPOD - Embedding Change Probabilities');

figure;
imshow(Stego-Cover,[]);
title('MiPOD - Changed Pixels (+1 -> white ,-1 -> black)');