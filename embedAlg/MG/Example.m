Cover = double(imread ('1.pgm'));

% Set the payload to 0.4 bpp
Payload = 0.4;

% MG embedding
tStart = tic;
[Stego, pChange, ChangeRate] = MG( Cover, Payload );
tEnd = toc(tStart);
fprintf('MG embedding is done in: %f (sec)\n',tEnd);

%%
figure;
imshow (Cover,[]);
title ('Cover image');

figure;
imshow(1-pChange/0.3333);
title('MG - Embedding Change Probabilities');

figure;
imshow(Stego-Cover,[]);
title('MG - Changed Pixels (+1 -> white ,-1 -> black)');