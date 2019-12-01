function ps_nr = cacul_psnr(img1, img2)
if(ischar(img1))
   img1 = imread(img1);
   img2 = imread(img2);
end
img1 = single(img1); img2 = single(img2);
[~, ps_nr] = psnr(img1(:,:,1),img2(:,:,1));
ps_nr = round(ps_nr,3);
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
