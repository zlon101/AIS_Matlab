function mextest()
% matlab coder mex function
cover = single(imread('E:\astego\Images\test\stegos\195.bmp'));
[stego, distortion] = HUGO(cover, single(0.4), single([1.1,1.2]));
