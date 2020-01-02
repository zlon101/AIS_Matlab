function result=MMD(FC,FS)
% 计算MMD
%%
t0=datetime('now');
% 标准化
[FC, mu, sigma] = zscore(FC, 1,1);
FS = bsxfun(@minus, FS, mu);
FS = bsxfun(@rdivide, FS, sigma+1e-30);
CMuSigma.mu = mu; CMuSigma.sigma = sigma;

old=''; N=size(FS,1); % 样本个数
CC = zeros(N); SS=zeros(N); CS=zeros(N);
for i=1:N
  Ci=FC(i,:); Si=FS(i,:);
	CC(i,:) = getK(Ci, FC);
  SS(i,:) = getK(Si, FS);
  CS(i,:) = getK(Ci, FS);
  % 打印
  msg=sprintf('- count: %3d/%d',i,N);
  fprintf([repmat('\b',1,length(old)),msg]); old=msg;
end
inds=eye(N,'logical');
CC(inds)=0; CS(inds)=0; SS(inds)=0;

%{
for i=1:N
  for j=1:N
    if(j==i) 
      continue;
    end
    Ci=FC(i,:); Cj=FC(j,:);
    Si=FS(i,:); Sj=FS(j,:);
    if(i>j)
      CC(i,j)= CC(j,i); SS(i,j)= SS(j,i);
    else
      CC(i,j)= getK(Ci,Cj);
      SS(i,j)= getK(Si,Sj);
    end
    CS(i,j)= getK(Ci,Sj);
  end
end
%}
s=sum(CC(:)) + sum(SS(:)) - 2*sum(CS(:));
result = sqrt( s / (N*(N-1)) );
fprintf('\n耗时：'); disp(datetime('now')-t0);
end

function r=getK(c,s)
  a= 1e-3;
  diff=s-c;
  r=zeros(1, size(diff,1));
  for i=1:size(diff,1)
    r(i) = exp(-a * norm(diff(i,:)).^2);
  end
end