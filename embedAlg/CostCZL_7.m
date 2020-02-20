function [rhoP1,rhoM1] = CostCZL7(cover)
% HUGO 代价函数
% 返回+1 -1 的代价
%% 
cover = single(cover);
wetCost = 10^8;
responseP1 = [-1, 1;];
% responseP1 = [0; 0; -1; +1; 0; 0];

% create mirror padded cover image
padSize = double(3);
coverPadded = padarray(cover, [padSize,padSize], 'symmetric');
% create residuals
CrezH = coverPadded(:, 1:end-1) - coverPadded(:, 2:end);

rhoM1 = zeros(size(cover),'single');  % declare cost of -1 change           
rhoP1 = zeros(size(cover),'single');  % declare cost of +1 change
t =1; a=2; % 滤波器阶数 & 权重
for row=1:size(cover, 1)
  for col=1:size(cover, 2)
    % Horizontal
    arrDH=[CrezH(row+2, col:col+t);
           CrezH(row+4, col:col+t)];
    rezH = CrezH(row+3, col:col+t);
    % CZL_7
    DCover = norm( reshape([arrDH;a*rezH],[],1) );
    rhoP1(row,col)= 1./(DCover+1e-20);
  end
end
% 平滑滤波
L= ones(9);
rhoP1= imfilter(rhoP1, L,'symmetric','conv','same')./sum(L(:));
rhoM1= rhoP1;

rhoM1(rhoM1>wetCost) = wetCost;
rhoP1(rhoP1>wetCost) = wetCost;
rhoP1(cover == 255) = wetCost;
rhoM1(cover == 0) = wetCost;
end
