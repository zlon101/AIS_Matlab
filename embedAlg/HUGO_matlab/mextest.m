function mextest()
% matlab coder mex function
root = 'E:\astego\Images\BOSS_ALL\';
name = '1013.pgm';
cover = single(imread([root,name]));
payload=single(0.2);
[rhoP1,rhoM1] = CostHUGO(cover);

% P = 1./rhoP1;
% P2 = filter2(ones(7), P,'same');
% nChange= round(numel(cover)*payload);
% clear payload
% HUGO_OPT_EMBEDS('E:\astego\Images\BOSS_ALL\', 'E:\astego\Images\HUGOOPT_04\',single(0.4));
% [stego, distortion] = HUGO(cover, payload);
% stego = EmbeddingSimulator(cover, rhoP1, rhoM1, message_length, false);