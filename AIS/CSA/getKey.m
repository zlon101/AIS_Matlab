function Key = getKey(Ab)
% Ab:  1*N array -- N个参数形式
% Key: char
%%
Key = "";
for i=1:length(Ab)
    Key = strjoin([Key,string(Ab(i))],'_');
end
Key = char(Key);
% Key = mat2str(Ab);
end