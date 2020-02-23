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
%{
rezH = cPadded(:, 1:end-1)- cPadded(:, 2:end);
rezV = cPadded(1:end-1, :)- cPadded(2:end, :);
rezD = cPadded(1:end-1, 1:end-1) - cPadded(2:end, 2:end);
rezMD= cPadded(1:end-1, 2:end) - cPadded(2:end, 1:end-1);
%}

T= 3; G=(T-1)*0.5;  % 滤波器阶数 & 权重, T=3,5,7
cH=1; cV=1;
rhoM1 = zeros(size(cover),'single');
rhoP1 = zeros(size(cover),'single');
for row=1:size(cover, 1)
  r=row+3;
  for col=1:size(cover, 2)
    c=col+3; % padSize=3;
    rs= r-G:r+G; cs= (c-G:c+G)-2; %-2
    
    subMatri= cPadded(rs,cs);
    resH= subMatri(:,1:end-1)-subMatri(:,2:end);
    resV= subMatri(1:end-1,:)-subMatri(2:end,:);
    
    resH(G+1,:)= resH(G+1,:).* 2;
    resV(:,G+1)= resV(:,G+1).* 2;
    resH=resH(:); resV=resV(:);
    
    rhoP1(row,col)= 1/(cH*norm(resH)+ cV*norm(resV)+ 1);
    
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
L= ones(9);
rhoP1= imfilter(rhoP1, L,'symmetric','conv','same')./sum(L(:));
% rhoP1 = ordfilt2(rhoP1,81,true(9),'symmetric');

rhoM1= rhoP1;
rhoM1(rhoM1>wetCost) = wetCost;
rhoP1(rhoP1>wetCost) = wetCost;
rhoP1(cover == 255) = wetCost;
rhoM1(cover == 0) = wetCost;
end