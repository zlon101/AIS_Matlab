function stego=detection(Tais,T,B)
% 用ensemble分类器提供的数据经过免疫算法运算,将处理后的特征用于检测,
% 验证免疫处理的效果

% -------------------------------------------------------------------------
% 修改C的值
%{
cover = load('cover.mat');
stego = load('stego.mat');
C=cover.F;

fs_ais=zeros(size(C));
nsamp=size(C,1);
for i=1:nsamp
    fs_ais(i,:)=C(i,:)*Tais*T+B;   %对应一幅载密图像特征
end
stego.F=fs_ais;
%}


cover = load('cover.mat');
stego = load('stego.mat');
names = intersect(cover.names,stego.names);             %交集
names = sort(names);

% Prepare cover features C,imgCove与imgSteg要对应
cover_names = cover.names(ismember(cover.names,names));
[cover_names,ix] = sort(cover_names);
C = cover.F(ismember(cover.names,names),:);
C = C(ix,:);

% Prepare stego features S
stego_names = stego.names(ismember(stego.names,names));
[stego_names,ix] = sort(stego_names);
S = stego.F(ismember(stego.names,names),:);
S = S(ix,:);

savaRelatCs(C,S,'E:\astego\算法对特征的影响\nsF5--CC-PEV\','CC-PEV');
% -------------------------------------------------------------------------
% 验证nsF5对载体-载密特征的影响
%{

ndw=size(C,2);              %特征维数
count=0;
while(count<=ndw)    
    close all;fh=figure;    
    set(fh,'visible','off');
    for i=1:4
        count=count+1;
        if(count>ndw)  
            break;
        end
        subplot(2,2,i);plot(C(:,count),S(:,count),'k.');
        title(['spam686-','多维特征点:',int2str(count)]);
        %hold on;plot(cf(:,i),polyval(p(i,:),cf(:,i)),'r-');    %对应拟合曲线            
    end    
    filename=[int2str(count),'.jpg'];
    saveas(fh,['.\feature\',filename]);                     %保存Figure 2窗口的图像
    if mod(count,10)==0
        fprintf('--%d:\n', count);
    end
end
%}

% -------------------------------------------------------------------------
% 验证免疫处理
fs_ais=zeros(size(C));
nsamp=size(C,1);
for i=1:nsamp
    fs_ais(i,:)=C(i,:)*Tais*T+B;   %对应一幅载密图像特征
end

% 计算fs与fc的相似度
d_c_s=zeros(nsamp,1);
for i=1:nsamp
    d_c_s(i)=pdist([S(i,:);C(i,:)],'euclidean');
end
d_c_s=mean(d_c_s);

% 计算fs_ais与fc的相似度
d_fsais=zeros(nsamp,1);
for i=1:nsamp
    d_fsais(i)=pdist([fs_ais(i,:);C(i,:)],'euclidean');
end
d_fs_ais=mean(d_fsais);

figure;
plot(S(:,1),'k.');hold on;plot(fs_ais(:,1),'bo');
legend({'原载密特征','免疫后载密特征'});

figure;
plot(C(:,1),'k.');hold on;plot(fs_ais(:,1),'bo');
legend({'载体特征','免疫后载密特征'});

fprintf('原始载体-载密特征距离:%.2f\n',d_c_s);
fprintf('免疫载体-载密特征距离:%.2f\n',d_fs_ais);
end
%}

%{
count=0;
TXT='';
while(count<4)
    TXT = updateTXT(TXT,sprintf(' - d_sub %d',count));
    count=count+1;
end  
%}