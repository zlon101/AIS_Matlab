function [cH,cV]= correl(x)
% 相关度: 相关性度量
resH= x(:,1:end-1)- x(:,2:end);
resV= x(1:end-1,:)- x(2:end,:);

stdH= std(resH(:));stdV= std(resV(:));
% varH= var(resH(:));varV= var(resV(:));
% nH= norm(resH(:));nV= norm(resV(:));
% mH= mean(resH(:));  mV= mean(resV(:));

cH= (stdH+stdV)/stdH;
cV= (stdH+stdV)/stdV;
cH=cH*0.5; cV=cV*0.5;
end

%             varH    varV    normH    normV    stdH    stdV    mH      mV
% 垂直条纹    419     0.6     1.3e3    49       20      0.7     0.78    -0.04                        
% 水平条纹    0.6     419     49       1.3e3    0.7     20      -0.04   0.78                     
% 
% 垂直1013    195     148     3.3e3    2.89e3   14      12      0.009   0.07
% 水平1013    152     232     1.6e3    2e3      12      15      0.009   0.07