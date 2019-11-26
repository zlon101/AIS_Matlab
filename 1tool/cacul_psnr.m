function ps_nr = cacul_psnr(imgName1, imgName2)
I1 = single(imread(imgName1));
I2 = single(imread(imgName2));
[~, ps_nr] = psnr(I1(:,:,1),I2(:,:,1));
end

%% 
function psnr=my_psnr(I0, I1)
[m,n]=size(I1);
mse=0;
[~,~,dim] = size(I0);
for i=1:dim
    I_0=I0(:,:,i);I_1=I1(:,:,i);
    ms=sum(sum( (I_0-I_1).^2 ));
    ms=ms/(m*n);
    mse=ms+mse;
end
mse=mse/dim;
psnr=10*log10( (255^2)/mse );
end
