function Key = caculMemoryKey(Ab)
% Ab:  1*N -- N个参数形式
%%
Key = "";
for i=1:length(Ab)
    Key = strjoin([Key,string(Ab(i))],'_');
end
% Key = mat2str(Ab);
end