function [stego, distortion] = HUGO(cover, payload, coefs)
% HUGO算法
% cover: single m*n array
% payload: single
%%
stego = zeros(size(cover),'single'); % 初始化 预分配内存
[rhoP1,rhoM1] = CostHUGO(cover, coefs);

%% Embedding
% embedding simulator - params.qarity \in {2,3}
stego = EmbeddingSimulator(cover, rhoP1, rhoM1, round(numel(cover)*payload), false);

%% compute distortion
distM1 = rhoM1(stego-cover==-1);
distP1 = rhoP1(stego-cover==1);
distortion = sum(distM1) + sum(distP1);
end
