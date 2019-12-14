function [fits,Mem]=fitOfAlg2(Abs,embedParas,Mem)
% second algorithm fit
% 失真函数的参数优化方案, 计算其适应度
%%
srcData = single(imread(embedParas.srcPath));
num = size(Abs,1);
fits = zeros(num,1,'single'); old='';
for i=1:num
  K = getKey(Abs(i,:));
  [~, ind] = ismember(K,Mem.K);
  if(ind>0)
    fits(i) = Mem.V(ind);
  else
    % 隐写
    stego = HUGO(srcData, embedParas.payLoad, Abs(i,:));
    fits(i) = calcuDist(srcData, stego);
    Mem.K{Mem.last}=K; Mem.V(Mem.last)=fits(i); Mem.last=Mem.last+1;
    
    if(Mem.last > length(Mem.V))
      last = Mem.last;
      VT=Mem.V;
      Mem.V = zeros(last+20,1,'single');
      Mem.V(1:last-1) = VT;
      Mem.last = last;
      clear VT;
    end
    % 打印
    msg=sprintf('- count: %3d/%d',i,num);
    fprintf([repmat('\b',1,length(old)),msg]);
    old=msg;
  end
% for-end  
end
clearvars -except fits Mem;
end