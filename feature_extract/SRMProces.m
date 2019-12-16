function F=SRMProces(S, reversal)
% reversal£ºÊÇ·ñ×ªÖÃ
propertys=fieldnames(S);
len=length(propertys);
F=[];
if(reversal)
  for i=1:len
    F=[F; getfield(S,propertys{i})];
  end
  F=F';
else
  for i=1:len
    F=[F, getfield(S,propertys{i})];
  end
end
end