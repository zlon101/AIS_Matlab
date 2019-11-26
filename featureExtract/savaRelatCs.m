function savaRelatCs(cf,sf,pathSavaImg,featureName,reduce)
% 保存载体-载密特征每一维的线性曲线
% cf,sf  特征向量
if(nargin<5)
    savaRelatCs_(cf,sf,pathSavaImg,featureName);
else
    savaRelatCs_error(cf,sf,pathSavaImg,featureName);
end    
end

function savaRelatCs_(cf,sf,pathSavaImg,featureName)
% 保存载体-载密特征每一维的线性曲线
% cf,sf  特征向量
ndw=size(cf,2);              %特征维数
count=0;

old='';
while(count<=ndw)    %ndw
    close all;fh=figure;    
    set(fh,'visible','off');
    for i=1:4
        count=count+1;
        if(count>ndw)  
            break;
        end
        subplot(2,2,i);plot(cf(:,count),sf(:,count),'k.');
        title([featureName,'-',int2str(count)]);
        xlabel(sprintf('Fc(:,%d)',count));
        ylabel(sprintf('Fs(:,%d)',count));
        %hold on;plot(cf(:,i),polyval(p(i,:),cf(:,i)),'r-');    %对应拟合曲线            
    end  
    filename=[int2str(count),'.jpg'];
    saveas(fh,[pathSavaImg,filename]);
    msg=sprintf('- count: %3d/%d',count,ndw);
    fprintf([repmat('\b',1,length(old)),msg]);
    old=msg;
end

%{
%-画出各维的拟合曲线
figure;
title('拟合曲线');xlabel('cover-img feature');ylabel('steg-img feature');
axis([0,1,0,1]);
for i=1:size(cf,2)
    hold on;
    plot(cf(:,i),polyval(p(i,:),cf(:,i)));
end
fprintf('曲线个数: %d', i);
%}
end

function savaRelatCs_error(cf,sf,pathSavaImg,featureName)
% 计算fs-fc，并保存

ndw=size(cf,2);              %特征维数
error=sf-cf;
count=0;

old='';
while(count<=ndw)    %ndw
    close all;fh=figure;    
    set(fh,'visible','off');
    for i=1:4
        count=count+1;
        if(count>ndw)  
            break;
        end
        subplot(2,2,i);plot(error(:,count),'k.');
        title([featureName,'-',int2str(count)]);
        xlabel(sprintf('Fc(:,%d)',count));
        ylabel(sprintf('Fs(:,%d)',count));
        %hold on;plot(cf(:,i),polyval(p(i,:),cf(:,i)),'r-');    %对应拟合曲线            
    end  
    filename=[int2str(count),'.jpg'];
    saveas(fh,[pathSavaImg,filename]);
    msg=sprintf('- count: %3d/%d',count,ndw);
    fprintf([repmat('\b',1,length(old)),msg]);
    old=msg;
end

end