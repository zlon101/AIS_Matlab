clear all;
global popsize length min max N code;
N=12;               % ÿ��Ⱦɫ�������ʮ���Ʊ���λ����
M=20;               % ��������
popsize=20;         %���ó�ʼ������Ⱥ���С
length=10;          % lengthΪÿ�λ���Ķ����Ʊ���λ��
chromlength=N*length;  %�ַ������ȣ����峤�ȣ���Ⱦɫ��Ķ����Ʊ��볤��
pc=0.7;                %���ý�����ʣ������н�������Ƕ�ֵ���������ñ仯�Ľ�����ʿ��ñ���ʽ��ʾ�����дһ��������ʺ�����������������ѵ���õ���ֵ��Ϊ�������
pm=0.3;               %���ñ�����ʣ�ͬ��Ҳ������Ϊ�仯��
bound={-100*ones(popsize,1),zeros(popsize,1)};min=bound{1};max=bound{2};
pop=initpop(popsize,chromlength);                     %���г�ʼ�����������������ʼȺ��
ymax=500;K=1;

%����C2:�������ͱ��룬ÿһ��Ϊһ�֣�code(1,:)��������code(2,:)��50����code(3,:)��150��
code =[-0.8180   -1.6201  -14.8590  -17.9706  -24.0737  -33.4498  -43.3949  -53.3849  -63.3451  -73.0295  -79.6806  -74.3230
       -0.7791   -1.2697  -14.8682  -26.2274  -30.2779  -39.4852  -49.4172  -59.4058  -69.3676  -79.0657  -85.8789  -81.0905
       -0.8571   -1.9871  -13.4385  -13.8463  -20.4918  -29.9230  -39.8724  -49.8629  -59.8215  -69.4926  -75.9868  -70.6706];
%ʵ�ʹ��ϲ�����ݱ��룬����Unnoralcode,188%
Unnoralcode=[-0.8864   -2.2743  -12.2676  -11.6813  -18.5298  -27.9828  -37.9340  -47.9246  -57.8820  -67.5433  -73.9248  -68.9759];
%%
for i=1:3   % 3�ֹ���ģʽ��ÿ��ģʽӦ�ò��� popsize �ּ���������壩��ÿ�ּ�����ĳ��Ⱥ͹��ϱ���ĳ�����ͬ
   for k=1:M
       [objvalue]=calobjvalue(pop,i);                 %����Ŀ�꺯��
       fitvalue=calfitvalue(objvalue); favg(k)=sum(fitvalue)/popsize;  %����Ⱥ����ÿ���������Ӧ��
       newpop=selection(pop,fitvalue); objvalue=calobjvalue(newpop,i); %ѡ��
       newpop=crossover(newpop,pc,k);  objvalue=calobjvalue(newpop,i); %����
       newpop=mutation(newpop,pm);     objvalue=calobjvalue(newpop,i); %����
       for j=1:N  %���룡
           temp(:,j)=decodechrom(newpop,1+(j-1)*length,length);      %��newpopÿ��(���壩ÿ�У�ÿ�λ���ת����ʮ������
           x(:,j)=temp(:,j)/(2^length-1)*(max(j)-min(j))+min(j);     % popsize��N ����ֵ���е���ת��Ϊ���������       
       end
       [bestindividual,bestfit]=best(newpop,fitvalue);%���Ⱥ������Ӧֵ���ĸ��弰����Ӧֵ
       if bestfit<ymax
          ymax=bestfit;K=k;
       end
       % y(k)=bestfit;
       if ymax<10     % ������ֵС���趨��ֵ��ֹͣ����
           X{i}=x;
           break
       end
       if k==1
           fitvalue_for=fitvalue;x_for=x;
       end
       result=resultselect(fitvalue_for,fitvalue,x_for,x);
       fitvalue_for=fitvalue;
       x_for=x;
       pop=newpop;
   end
   X{i}=result;   % ��i����ϵ�popsize�������
   distance=0;    % ���濪ʼ���� Unnoralcode ����ÿһ����ϵĸ��ʣ�����Ӧ��
   for j=1:N
       distance=distance+(result(:,j)-Unnoralcode(j)).^2;   % ����popsize����ͬ�ľ���
   end
   distance=sqrt(distance);D=0;
   for p=1:popsize
       if distance(p)<80 % Ԥ�跧ֵ
           D=D+1;
       end
   end
   P(i)=D/popsize; %Unnoralcode����ÿ�ֹ������͵ĸ��ʣ������ռ������-->P!     
end

X;   % ���Ϊ(i*popsie)������������壩
plot(1:M,favg)