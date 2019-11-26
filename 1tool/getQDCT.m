function [QDCT, quality, QTable, Cb, Cr]=getQDCT(I, quality)
% I:像素矩阵或图像文件名,若是文件名，则调用CCPEV中的函数
% 计算QDCT系数
% quality=50时，量化表为标准量化表
% quality=100时，量化表为1，没有量化
if(ischar(I))
    [jobj, QTable, quality] = DCTPlane(I);
    Cb=0;Cr=0;
    QDCT = jobj;
    % QDCT = jobj.coef_arrays{1};
else
    I=double(I)-128;
    if(length(size(I)) > 2)
        R=I(:,:,1);G=I(:,:,2);B=I(:,:,3);
        [Y,Cb,Cr] = rgb2yCbCr( R,G,B );
        I = Y;
        %disp('--彩色图像--');
    else
        Cb=0;Cr=0;
    end

    QTable=calcuQTable(quality);
    % DCT变换
    dct_fun=@(block_struct) dct2(block_struct.data, [8,8]);
    coef_dcted=(blockproc(I, [8,8], dct_fun));
    % 量化 ok
    quan_fun=@(block_struct) ( block_struct.data ./ QTable); 
    QDCT = round(blockproc(coef_dcted, [8 8], quan_fun));
end
end

% 计算对应量化质量的量化表
function QTable=calcuQTable(quality)
% 
QTable=[16 11 10 16 24 40 51 61                 %Y分量量化系数
    12 12 14 19 26 58 60 55
    14 13 16 24 40 57 69 56
    14 17 22 29 51 87 80 62
    18 22 37 56 68 109 103 77
    24 35 55 64 81 104 113 92
    49 64 78 87 103 121 120 101
    72 92 95 98 112 100 103 99];

QTable=double(QTable);quality=double(quality);
if(0<quality && quality<50)
    for i=1:size(QTable,1)
        for j=1:size(QTable,2)
            QTable(i,j)=QTable(i,j)*50/quality;
        end
    end
else
    for i=1:size(QTable,1)
        for j=1:size(QTable,2)
            QTable(i,j)=max( QTable(i,j)*(2-quality/50), 1);
        end
    end
end
QTable = round(QTable);
end

% jpeg_read提取QDCT系数
function [jobj, q_table, quality] = DCTPlane(path)
    jobj=jpeg_read(path);
    QDCT=jobj.coef_arrays{1};
    q_table=jobj.quant_tables{1};
    
    % 计算Quality
    Q100 = q_table(end,6);
    if(Q100==1)
        quality = 100;
    elseif(Q100 < 100)
        quality = ( 2 - (Q100*0.01) ) * 50;
    else
        quality = (50*100)/Q100;
    end    
end