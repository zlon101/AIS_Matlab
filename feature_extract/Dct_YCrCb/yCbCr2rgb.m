function [ R,G,B ] = yCbCr2rgb( Y,Cb,Cr )
%UNTITLED2 此处显示有关此函数的摘要

R=Y+1.402*Cr ;
G=Y-0.34414*Cb -0.71414*Cr;
B=Y +1.772*Cb;             

end

