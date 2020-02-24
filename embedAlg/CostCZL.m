function [rhoP1,rhoM1] = CostCZL(cover, TFilter)
% HUGO 代价函数
% 返回+1 -1 的代价
%% 
cover = single(cover);
wetCost = 10^8;
% create mirror padded cover image
padSize = double(3);
cPadded = padarray(cover, [padSize,padSize], 'symmetric');
% create residuals
%{
rezH = cPadded(:, 1:end-1)- cPadded(:, 2:end);
rezV = cPadded(1:end-1, :)- cPadded(2:end, :);
rezD = cPadded(1:end-1, 1:end-1) - cPadded(2:end, 2:end);
rezMD= cPadded(1:end-1, 2:end) - cPadded(2:end, 1:end-1);
%}

%% 分割,计算不同方向的相关度
gap= 128;
nBlockH= size(cover,1)/gap; nBlockV=size(cover,2)/gap;
rows= 0:gap:size(cover,1);
cols= 0:gap:size(cover,2);
cH= zeros(nBlockH,nBlockV); cV=zeros(nBlockH,nBlockV);
for i=1:nBlockH
  r= rows(i)+1 : rows(i+1);
  for j=1:nBlockV
    c= cols(j)+1 : cols(j+1);
    subImg= cover(r, c);
    [cH(i,j),cV(i,j)]= correl(subImg);
  end
end

%% distortion cost
T= 3; G=(T-1)*0.5;  % T阶领域, T=3,5,7
% cH=1; cV=1;
rhoM1 = zeros(size(cover),'single');
rhoP1 = zeros(size(cover),'single');
for row=1:size(cover, 1)
  i=ceil(row/gap);
  r= row+3;
  for col=1:size(cover, 2)
    j= ceil(col/gap);
    c=col+3; % padSize=3;
    rs= r-G:r+G; cs= (c-G:c+G);% -2; %偏移
    
    subMatri= cPadded(rs,cs);
    resH= subMatri(:,1:end-1)-subMatri(:,2:end);
    resV= subMatri(1:end-1,:)-subMatri(2:end,:);
    
    resH(G+1,:)= resH(G+1,:).* 2;
    resV(:,G+1)= resV(:,G+1).* 2;
    resH=resH(:); resV=resV(:);
    
    rhoP1(row,col)= 1/(cH(i,j)*norm(resH)+ cV(i,j)*norm(resV)+ 1);
    
    % CZL8 T=3
    %{
    subH= [rezH(r-1, c-T:c+T-1); 2.* rezH(r, c-T:c+T-1); rezH(r+1, c-T:c+T-1)];
    subV= [rezV(r-T:r+T-1, c-1); 2.* rezV(r-T:r+T-1, c); rezV(r-T:r+T-1, c+1)];
    rhoP1(row,col)= 1/ (cH*norm(subH(:))+ cV*norm(subV(:))+ 1e-20);
    %}
    % CZL7
    %{
    %DCover = norm( reshape([arrDH;a*rezH],[],1) );
    %rhoP1(row,col)= 1./(DCover+1e-20);
    %}
  end
end

% 平滑滤波
% TFilter = 15;
L= ones(TFilter);
rhoP1= imfilter(rhoP1, L,'symmetric','conv','same')./sum(L(:));
% rhoP1 = ordfilt2(rhoP1,81,true(9),'symmetric');

rhoM1= rhoP1;
rhoM1(rhoM1>wetCost) = wetCost;
rhoP1(rhoP1>wetCost) = wetCost;
rhoP1(cover == 255) = wetCost;
rhoM1(cover == 0) = wetCost;
end