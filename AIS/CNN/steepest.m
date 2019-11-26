function [k ender]=steepest(f,x,e)
%梯度下降法,f为目标函数（两变量x1和x2），x为初始点,如[3;4]
syms x1 x2 m; %m为学习率
d=-[diff(f,x1);diff(f,x2)];  %分别求x1和x2的偏导数，即下降的方向
flag=3;  %循环标志
k=0; %迭代次数
while(flag)
    d_temp=subs(d,x1,x(1));      %将起始点代入，求得当次下降x1梯度值
    d_temp=subs(d_temp,x2,x(2)); %将起始点代入，求得当次下降x2梯度值
    nor=norm(d_temp); %范数
    if(nor>=e)
        x_temp=x+m*d_temp;            %改变初始点x的值
        f_temp=subs(f,x1,x_temp(1));  %将改变后的x1和x2代入目标函数
        f_temp=subs(f_temp,x2,x_temp(2));
        h=diff(f_temp,m);  %对m求导，找出最佳学习率
        m_temp=solve(h);   %求方程，得到当次m
        x=x+m_temp*d_temp; %更新起始点x
        k=k+1;
    else
        flag=0;
    end
    flag=flag-1;
end
ender=double(x);  %终点
end