function [stego, distortion] = HUGO_Input(cover, payload)
% HUGO算法，入口文件
% author: czl
% 
%% 
% set params
params.gamma = 1;
params.sigma = 1;

[stego, distortion] = HUGO_like(cover, payload, params);
end