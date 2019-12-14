% function mextest()
% matlab coder mex function
cover = single(imread('E:\astego\Images\Experis\195.pgm'));
payload=single(0.2);
message_length = round(numel(cover)*payload);

HUGO_OPT_EMBEDS('E:\astego\Images\BOSS_ALL\', 'E:\astego\Images\HUGOOPT_04\',single(0.4));

% [stego, distortion] = HUGO(cover, payload);
% imshow(stego-cover,[]);
% [rhoP1,rhoM1] = CostHUGO(cover);
% figure('name','HUGO_+1¸ÅÂÊ'); imshow(rhoP1.\1,[]);

% stego = EmbeddingSimulator(cover, rhoP1, rhoM1, message_length, false);
% [stego, distortion] = HUGO(cover, single(0.4), single([1.1,1.2]));

