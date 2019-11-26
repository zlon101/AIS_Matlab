function [fsByT,FsByAis]=verifyAis(fc,fs,Tais,T,B,isSava,featname,pathSavaImg)
% 测试T和B的正确性,计算FsAis=Fc*Tais*T
% 验证Tais的效果
if (nargin<6)
    isSava = 0;
end

setSampNum=200;                         % 设置显示的样本数
nsamp=size(fc,1);

if(setSampNum>nsamp)
    setSampNum=nsamp;
end
% setSampNum=nsamp;

ndw=size(fc,2);
fsByT=zeros(nsamp,ndw);
for i=1:nsamp
    fsByT(i,:)=fc(i,:)*T+B;            % Fc*T+B=Fs,一行代表一个样本观测值
end
[coeef,score0,latent,~]=pca(fs);   %pca降维
pcaFs=score0(:,1:2);
[coeef,score,latent,~]=pca(fsByT);
pcaFsByT=score(:,1:2);
% 计算fs与fsByT的相似度
D1=zeros(size(fs,1),1);
for i=1:size(fs,1)
    D1(i)=pdist([fsByT(i,:);fs(i,:)],'euclidean');
end
dFs_FsT=mean(D1);

% 保存各个维度下,Te的效果
if(isSava)
old='';
count=0;                            % 第count维度分量
while(count<=ndw)
    close all;fh=figure;
    set(fh,'visible','off');    
    count=count+1;
    if(count>ndw)
        break;
    end        
    plot(fs(1:setSampNum,count),'k.');hold on;plot(fsByT(1:setSampNum,count),'bo');        
    legend({'Fs','Fc*T'});
    title('Fc*T与Fs的相似度');
    
    filename=[featname,'-',int2str(count),'.jpg'];
    saveas(fh,[pathSavaImg,filename]);
    msg=sprintf('- count: %3d/%d',count,ndw);
    fprintf([repmat('\b',1,length(old)),msg]);
    old=msg;
end
end
% saveas(f1,'转换矩阵对隐写算法的模拟效果_spam.jpg');


%% 比较免疫处理后的特征FsByAis与 Fc相似度
FsByAis=zeros(nsamp,ndw);
for i=1:nsamp
    FsByAis(i,:)=fc(i,:)*Tais*T+B;   %对应一幅载密图像特征
end
% 计算fs_ais与fc的相似度
D2=zeros(size(fs,1),1);
for i=1:size(fs,1)
    D2(i)=pdist([FsByAis(i,:);fc(i,:)],'euclidean');
end
dFs_Ais=mean(D2);
fprintf('distance of Fc*T and  Fs: %d\ndistance of Fc*Tais*T and Fc: %d\n', dFs_FsT,dFs_Ais);

% Fc*Tais*T与Fc的相似度
figure; % suptitle('Fc*Tais*T与Fc的相似度')
subplot(1,2,1); 
plot(fs(1:setSampNum,1),'k.');hold on;
plot(fsByT(1:setSampNum,1),'bo');
title('Fc*T & Fs in dimension 1');legend({'Fs','Fc*T'});
subplot(1,2,2); 
plot(fs(1:setSampNum,2),'k.');hold on;
plot(fsByT(1:setSampNum,2),'bo');
title('Fc*T & Fs in dimension 2');legend({'Fs','Fc*T'});

%% 保存各个维度下,免疫的效果
if(0)
count=0;
choice_dw=randi([1 ndw],4,1);           %4*1 修改
while(count<=ndw)
    close all;fh=figure;    
    set(fh,'visible','off');
    for i=1:4
        count=count+1;
        if(count>ndw)
            break;
        end
        subplot(2,2,i);
        %plot(fc(1:setSampNum,count),'k+');hold on;plot(FsByAis(1:setSampNum,count),'bo');
        plot(fc(1:setSampNum,choice_dw(count)),'k.');hold on;plot(FsByAis(1:setSampNum,choice_dw(count)),'bo');
        legend({'Fc','Fc*Tais*T'});
        title([featname,'--:',int2str(count)]);        
    end    
    filename=[featname,'-',int2str(count),'.jpg'];
    saveas(fh,[pathSavaImg,filename]);
    break;                          % 随机保存其中4个维度
end
end