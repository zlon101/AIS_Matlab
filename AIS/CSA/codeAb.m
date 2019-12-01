function gene = codeAb(Ab,vmin,len,precis)
% 输入抗体，输出基因序列
for i=1:length(Ab)
  str = fliplr(dec2bin( (Ab(i)-vmin)/precis ));
  if(length(str)<len)
    s = mat2str(zeros(1,len-length(str)));
    str = [str,s];
  end
end
gene = str(str~=' ' & str~='[' & str~=']');
end

