cRoot= 'E:\astego\Images\3\cover\';
format= 'bmp';
name= 'mandril_gray'; name= [name,'.',format];
cpath = [cRoot,name];
payload=single(0.4);
cover= single(imread(cpath));

sharpRoot= 'E:\astego\Images\3\sharpImg\';
stegoRoot= 'E:\astego\Images\3\sharpStego\';
sharpImg= sharpen(cover, 0.7);
stego= HUGO_like(uint8(sharpImg), payload);
imwrite(uint8(sharpImg),[sharpRoot,name],format);
imwrite(uint8(stego),[stegoRoot,name],format);

% [rhoP1,rhoM1] = CostUNWD(cover);
% [rhoP1,rhoM1] = CostCZL(cover);
% stego= embedAlgCZL(cpath,payload,'1');

%---------------------------------------------------------------------------
% rhoP1= rhoP1(154:291,58:228);
% P=1./rhoP1;
% figure('name','HUGO');histogram(P1);
% figure('name','HUGO');imshow(P1,'DisplayRange',[0,1.5],'Border','tight');
% figure('name','CZL');histogram(P2);
% figure('name','CZL');imshow(P2,'DisplayRange',[1,50],'Border','tight');
%---------------------------------------------------------------------------
% sigma= 13;
% Fsize= 2*ceil(2*sigma)+1;
% Fsize=13;
% Filter1 = fspecial('gaussian',Fsize,sigma);
%---------------------------------------------------------------------------

%---------------------------------------------------------------------------