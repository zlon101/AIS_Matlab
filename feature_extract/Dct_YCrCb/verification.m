function coef_quaned=verification(I,Q)
% 验证rgb2YcrCb
% 验证DCT

% 标准量化表
Qy=[16 11 10 16 24 40 51 61                 %Y分量量化系数
    12 12 14 19 26 58 60 55
    14 13 16 24 40 57 69 56
    14 17 22 29 51 87 80 62
    18 22 37 56 68 109 103 77
    24 35 55 64 81 104 113 92
    49 64 78 87 103 121 120 101
    72 92 95 98 112 110 103 99];
pixel_test=[52 55 61 66  70  61 64 73
            63 59 55 90  109 85 69 72
            62 59 68 113 144 104 66 73
            63 58 71 122 154 106 70 69
            67 61 68 104 126 88 68 70
            79 65 60 70 77 68 58 75
            85 71 64 59 55 61 65 83
            87 79 69 68 65 76 78 94];
rank=80;
if (nargin<2)
    Q_ = Qy*100/rank;
else
    Q_=Q{1};
end
% Q_=ones(8,8);
%pixel_test=pixel_test-128;

imagname='D:\MATLAB_Software\myInstall\bin\image\lena_color_512.jpg';
%I=double(imread(imagname));
% I=I(1:64,1:64,:);
if(length(size(I))>2)               %彩色图像
    R=I(:,:,1);G=I(:,:,2);B=I(:,:,3);
    R=R-128;G=G-128;B=B-128;
    [Y,Cb,Cr]=rgb2yCbCr(R,G,B);
else
    Y=I;
end
% DCT变换
dct_fun=@(block_struct)dct2(block_struct.data,[8,8]);        %函数
coef_dcted=(blockproc(Y,[8 8],dct_fun));

% 量化 ok
quan_fun=@(block_struct) ( block_struct.data./Q_); 
coef_quaned=blockproc(coef_dcted,[8 8],quan_fun);           
% statistic(coef_quaned(:),imagname);

%%
while(0)
% ais处理
coef_aised=ais(round(coef_quaned));
% coef_aised=(coef_quaned);
% statistic(round(coef_quaned(:)),round(coef_aised(:)),'**');


% 验证空间转换和DCT是否正确
iquan_fun=@(block_struct) ( block_struct.data.*Q_); 
coef_quaned_recon=blockproc(coef_aised,[8 8],iquan_fun);
% DCT反变换
idct_fun=@(block_struct)idct2(block_struct.data,[8,8]);
Y_2=(blockproc(coef_quaned_recon,[8 8],idct_fun)); 

[ R2,G2,B2 ] = yCbCr2rgb( Y_2,Cb,Cr );
R2=R2+128;G2=G2+128;B2=B2+128;
I2(:,:,1)=double(R2);
I2(:,:,2)=double(G2);
I2(:,:,3)=double(B2);
imshow(uint8(I2));

% 提取特征
F=spam_extract_2(I2,3);
% if()
% end
end
end

function F = spam_extract_2(X,T)

% horizontal left-right
D = X(:,1:end-1) - X(:,2:end);
L = D(:,3:end); C = D(:,2:end-1); R = D(:,1:end-2);
Mh1 = GetM3(L,C,R,T);

% horizontal right-left
D = -D;
L = D(:,1:end-2); C = D(:,2:end-1); R = D(:,3:end);
Mh2 = GetM3(L,C,R,T);
% vertical bottom top
D = X(1:end-1,:) - X(2:end,:);
L = D(3:end,:); C = D(2:end-1,:); R = D(1:end-2,:);
Mv1 = GetM3(L,C,R,T);

% vertical top bottom
D = -D;
L = D(1:end-2,:); C = D(2:end-1,:); R = D(3:end,:);
Mv2 = GetM3(L,C,R,T);

% diagonal left-right
D = X(1:end-1,1:end-1) - X(2:end,2:end);
L = D(3:end,3:end); C = D(2:end-1,2:end-1); R = D(1:end-2,1:end-2);
Md1 = GetM3(L,C,R,T);

% diagonal right-left
D = -D;
L = D(1:end-2,1:end-2); C = D(2:end-1,2:end-1); R = D(3:end,3:end);
Md2 = GetM3(L,C,R,T);

% minor diagonal left-right
D = X(2:end,1:end-1) - X(1:end-1,2:end);
L = D(1:end-2,3:end); C = D(2:end-1,2:end-1); R = D(3:end,1:end-2);
Mm1 = GetM3(L,C,R,T);

% minor diagonal right-left
D = -D;
L = D(3:end,1:end-2); C = D(2:end-1,2:end-1); R = D(1:end-2,3:end);
Mm2 = GetM3(L,C,R,T);

F1 = (Mh1+Mh2+Mv1+Mv2)/4;
F2 = (Md1+Md2+Mm1+Mm2)/4;
F = [F1;F2];
end

function M = GetM3(L,C,R,T)
% marginalization into borders
L = L(:); L(L<-T) = -T; L(L>T) = T;
C = C(:); C(C<-T) = -T; C(C>T) = T;
R = R(:); R(R<-T) = -T; R(R>T) = T;

% get cooccurences [-T...T]
M = zeros(2*T+1,2*T+1,2*T+1);
for i=-T:T
    C2 = C(L==i);
    R2 = R(L==i);
    for j=-T:T
        R3 = R2(C2==j);
        for k=-T:T
            M(i+T+1,j+T+1,k+T+1) = sum(R3==k);
        end
    end
end

% normalization
M = M(:)/sum(M(:));
end