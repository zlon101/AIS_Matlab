function [rhoP1,rhoM1] = CostCZL(cover)
% HUGO 代价函数
% 返回+1 -1 的代价
%% 
cover = single(cover);
wetCost = 10^8;

% create mirror padded cover image
padSize = double(3);
cPadded = padarray(cover, [padSize,padSize], 'symmetric');
% create residuals
rezH = cPadded(:, 1:end-1)- cPadded(:, 2:end);
rezV = cPadded(1:end-1, :)- cPadded(2:end, :);
rezD = cPadded(1:end-1, 1:end-1) - cPadded(2:end, 2:end);
rezMD= cPadded(1:end-1, 2:end) - cPadded(2:end, 1:end-1);

rhoM1 = zeros(size(cover),'single');  % declare cost of -1 change           
rhoP1 = zeros(size(cover),'single');  % declare cost of +1 change
T= 5; T=(T-1)*0.5;  % 滤波器阶数 & 权重, T=3,5,7
cH=1; cV=1;
for row=1:size(cover, 1)
  r=row+3;
  for col=1:size(cover, 2)
    c=col+3; % padSize=3;
    % CZL09 T=5
    %
    subH= 2.* rezH(r, c-T:c+T-1);
    subV= 2.* rezV(r-T:r+T-1, c);
    subD= [rezD(r-2,c-2); rezD(r-1,c-1); rezD(r,c); rezD(r+1,c+1);];
    subMD=[rezMD(r-2,c+1); rezMD(r-1,c); rezMD(r,c-1); rezMD(r+1,c-2);];
    rhoP1(row,col)= 1/(norm(subH(:))+ norm(subV(:))+ norm(subD(:))+ norm(subMD(:))+ 1e-20);
    %}
    
    % CZL8 T=3
    %{
    subH= [rezH(r-1, c-T:c+T-1); 2.* rezH(r, c-T:c+T-1); rezH(r+1, c-T:c+T-1)];
    subV= [rezV(r-T:r+T-1, c-1); 2.* rezV(r-T:r+T-1, c); rezV(r-T:r+T-1, c+1)];
    rhoP1(row,col)= 1/ (cH*norm(subH(:))+ cV*norm(subV(:))+ 1e-20);
    %}
    
    % CZL7
    %DCover = norm( reshape([arrDH;a*rezH],[],1) );
    %rhoP1(row,col)= 1./(DCover+1e-20);
  end
end

% 平滑滤波
L= ones(9);
rhoP1= imfilter(rhoP1, L,'symmetric','conv','same')./sum(L(:));
% rhoP1 = ordfilt2(rhoP1,81,true(9),'symmetric');

rhoM1= rhoP1;
rhoM1(rhoM1>wetCost) = wetCost;
rhoP1(rhoP1>wetCost) = wetCost;
rhoP1(cover == 255) = wetCost;
rhoM1(cover == 0) = wetCost;
end