function [Y,Cb,Cr] = rgb2yCbCr(R,G,B)
%rgb×ªµ½YCrCb

Y=0.299*R+0.587*G+0.114*B;
Cb=-0.1687*R-0.3313*G+0.5*B; 
Cr=0.5*R-0.4187*G-0.0813*B;
end

