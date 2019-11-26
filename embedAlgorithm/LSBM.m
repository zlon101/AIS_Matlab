function steg =LSBM(I, payLoad)
% LSB Match算法
% I:像素矩阵

[R,C]=size(I);
I=double(I);
% 置乱I
% rng(7);
randp=randperm(R*C);
randp=reshape(randp, [R,C]);
steg=I(randp);
Ibak=steg;
msg=double( rand(round(payLoad*R*C), 1)<0.5);
pMsg=1;
isEnd=0;
for i=1:R
    if(isEnd || i+1>R)
        break;
    end
    
    for j=1:2:C
        if(j+1>C)
            break;
        end
        x1=steg(i,j);  x2=steg(i,j+1);
        if(x1==255 || x1==0 || x2==255 || x2==0)
            continue;
        end
        
        if( LSB(x1) == LSB( msg(pMsg)) )
            if( LSB(msg(pMsg+1)) ~= subfunc(x1, x2) )
                if( msg(pMsg+1) )
                    steg(i,j+1)=steg(i,j+1)+1;
                else
                    steg(i,j+1)=steg(i,j+1)-1;
                end
            end                    
        else
            if( LSB(msg(pMsg+1)) == subfunc(x1-1,x2))
                steg(i,j)=steg(i,j)-1;
            else
                steg(i,j)=steg(i,j)+1;
            end
        end
        pMsg=pMsg+2;
        if(pMsg > length(msg)-1)
            % fprintf('\n---embedding end---\n');
            isEnd=1;
            break;
        end
    end
end
%fprintf('%d\n', pMsg);
% 提取
% msgEx=extract(steg, length(msg));
% err=msg-msgEx;

steg(randp)=steg;
steg=round(reshape(steg,[R,C]));
steg = double(uint8(steg));
end

function result=subfunc(x1, x2)
% 子函数,匹配函数
% 
x1=double(x1);x2=double(x2);
result=floor(x1/2)+x2;
result=LSB(result);
end

function LSB=LSB(D)
% 计算LSB
LSB=mod(D,2);
end

function msg=extract(I, len)
% 提取
% 
[R,C]=size(I);
msg=zeros(len,1);
pmsg=1;
isEnd=0;
for i=1:R
    if(isEnd || i+1>R)
        break;
    end
    
    for j=1:2:C
        if(isEnd || j+1>C)
            break;
        end
        x1=I(i,j);  x2=I(i,j+1);
%         if(x1==255 || x1==0 || x2==255 || x2==0)
%             continue;
%         end
        
        m1=LSB(x1);        
        m2=subfunc(x1,x2);
        msg(pmsg)=m1;
        msg(pmsg+1)=m2;
        pmsg=pmsg+2;
        if(pmsg>len)
            isEnd=1;
        end
    end
end
fprintf('\n---extract end---\n');
end