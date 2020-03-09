function [rhoP1,rhoM1] = CostCZL(cover)
% HUGO 代价函数
% 返回+1 -1 的代价
%% 
% disp('CZL7!');
wetCost = 10^8;
cover = single(cover);
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
%{
gap= 64;
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
%}
%% distortion cost
a= 1;          %cH=1; cV=1; 
T= 3; G=(T-1)*0.5;  % T阶领域, T=3,5,7
rhoM1 = zeros(size(cover),'single');
rhoP1 = zeros(size(cover),'single');
% optP1= zeros(size(cover),'logical'); optM1= zeros(size(cover),'logical');
for row=1:size(cover, 1)
  r= row+3;
  for col=1:size(cover, 2)
    c=col+3;
    rs= r-G:r+G; cs= (c-G:c+G); % -2; %偏移
    subMatri= cPadded(rs,cs);
    % 残差
    resH= subMatri(:,1:end-1)-subMatri(:,2:end);
    resV= subMatri(1:end-1,:)-subMatri(2:end,:);
    %resD= subMatri(1:end-1, 1:end-1)- subMatri(2:end, 2:end);
    %resMD=subMatri(1:end-1, 2:end)- subMatri(2:end, 1:end-1);

    resH(G+1,:)= resH(G+1,:).* a;
    resV(:,G+1)= resV(:,G+1).* a;
    resH=resH(:); resV=resV(:); % resD=resD(:); resMD=resMD(:);
    %resH=resH(G+1,:); resV=resV(:,G+1);

    %增加对角
    %{
    resMD=[subMatri(1,1)-subMatri(2,2);
           subMatri(2,2)-subMatri(3,3);
           subMatri(3,3)-subMatri(4,4);
           subMatri(4,4)-subMatri(5,5);];
    resD= [subMatri(1,5)-subMatri(2,4);
           subMatri(2,4)-subMatri(3,3);
           subMatri(3,3)-subMatri(4,2);
           subMatri(4,2)-subMatri(5,1);];
    resD=resD(:); resMD=resMD(:);
    %rhoP1(row,col)= 1/(min([norm(resH),norm(resV),norm(resMD),norm(resD)]) +  1);
    %}
    x= [norm(resH); norm(resV)]; % ; norm(resD); norm(resMD)
    rhoP1(row,col)= 1/(min(x)+ 1);
    %}
    
    %非对称嵌入; rhoP1!=rhoM1;
    %{
    resP1H= resHmatri(subP1); resP1V= resVmatri(subP1);
    resM1H= resHmatri(subM1); resM1V= resVmatri(subM1);
    resP1H=resP1H(:); resP1V=resP1V(:); resM1H=resM1H(:); resM1V=resM1V(:);
    
    rhoP1(row,col)= 1/(cH(i,j)*norm(resP1H) + cV(i,j)*norm(resP1V)+ 1);
    rhoM1(row,col)= 1/(cH(i,j)*norm(resM1H) + cV(i,j)*norm(resM1V)+ 1);
    %}
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

%% 平滑滤波
TFilter = 11;
L= ones(TFilter);  %L = fspecial('gaussian',TFilter,sigma);
rhoP1= imfilter(rhoP1, L,'symmetric','conv','same')./sum(L(:));

% sigma = 13; % sigma越大,越接近均值 FSize=2*ceil(2*sigma)+1
% rhoP1=imgaussfilt(rhoP1,sigma,'Padding','symmetric');
% rhoP1 = ordfilt2(rhoP1,TFilter^2,true(TFilter),'symmetric');

rhoM1= rhoP1;
rhoM1(rhoM1>wetCost) = wetCost;
rhoP1(rhoP1>wetCost) = wetCost;
rhoP1(cover == 255) = wetCost;
rhoM1(cover == 0) = wetCost;

%{
  function resH=resHmatri(x)
    resH= x(:,1:end-1)- x(:,2:end);
    resH(G+1,:)= resH(G+1,:).* 2;
  end
  function resV=resVmatri(x)
    resV= x(1:end-1,:)- x(2:end,:);
    resV(:,G+1)= resV(:,G+1).* 2;
  end
%}
end