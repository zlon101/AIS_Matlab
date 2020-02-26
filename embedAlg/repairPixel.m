function [optP1,optM1] = repairPixel(cover)
% 找到具有弥补相关性的像素
%% 
cover = single(cover);
padSize = double(3);
cPadded = padarray(cover, [padSize,padSize], 'symmetric');

T= 3; G=(T-1)*0.5;  % T阶领域, T=3,5,7
optP1= zeros(size(cover),'logical');
optM1= zeros(size(cover),'logical');
for row=1:size(cover, 1)
  r= row+3;
  for col=1:size(cover, 2)
    c=col+3; % padSize=3;
    rs= r-G:r+G; cs= (c-G:c+G);
    
    subMatri= cPadded(rs,cs);
    %{
    vcenter=subMatri(G+1,G+1); subMatri(G+1,G+1)=nan;
    x=sort(subMatri(:));
    if( var(x(1:end-1))==0 && vcenter~=x(1))
      optP1(row,col)=1;
    end
    %}
   
    % 弥补嵌入, 嵌入后的相关度更高
    subP1= subMatri; subM1= subMatri;
    subP1(G+1,G+1)= subP1(G+1,G+1)+1;
    subM1(G+1,G+1)= subM1(G+1,G+1)-1;
    subP1=subP1(:); subM1=subM1(:);
    
    if(var(subP1)==0)
      optP1(row,col)=1;
    end
    if(var(subM1)==0)
      optM1(row,col)=1;
    end
    %}
  end
end