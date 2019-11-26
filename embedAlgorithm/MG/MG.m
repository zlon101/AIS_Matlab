function [Stego, pChange, ChangeRate] = MG ( Cover, Payload )
% -------------------------------------------------------------------------
% Multivariate Gaussian Embedding    |   September 2015    |   version 1.0 
% -------------------------------------------------------------------------
% INPUT:
%  - Cover - Path to the cover image or the cover image itself.
%  - Payload - Embedding payload in bits per pixel (bpp).
% OUTPUT:
%  - Stego - Resulting image with embedded payload
%  - pChange - Embedding change probabilities. 
%  - ChangeRate - Average number of changed pixels
% -------------------------------------------------------------------------
% Copyright (c) 2015 DDE Lab, Binghamton University, NY.
% All Rights Reserved.
% -------------------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software for
% educational, research and non-profit purposes, without fee, and without a
% written agreement is hereby granted, provided that this copyright notice
% appears in all copies. The program is supplied "as is," without any
% accompanying services from DDE Lab. DDE Lab does not warrant the
% operation of the program will be uninterrupted or error-free. The
% end-user understands that the program was developed for research purposes
% and is advised not to rely exclusively on the program for any reason. In
% no event shall Binghamton University or DDE Lab be liable to any party
% for direct, indirect, special, incidental, or consequential damages,
% including lost profits, arising out of the use of this software. DDE Lab
% disclaims any warranties, and has no obligations to provide maintenance,
% support, updates, enhancements or modifications.
% -------------------------------------------------------------------------
% Contact: vsedigh1@binghamton.edu | fridrich@binghamton.edu
%          September 2015
%          http://dde.binghamton.edu/download/
% -------------------------------------------------------------------------
% References:
% [1] - J. Fridrich and J. Kodovsky. Multivariate Gaussian model for 
% designing additive distortion for steganography. Proc. IEEE, ICASSP, 
% Vancouver, Canada, May 26-31, 2013.
% -------------------------------------------------------------------------

% Read and convert the input cover image into double format
if ischar( Cover )
    Cover = double( imread(Cover) );
else
    Cover = double( Cover );
end

% Compute Variance and do the flooring for numerical stability
Variance = VarianceEstimation(Cover);
Variance(Variance< 1) = 1;

% Compute Fisher information and smooth it
FisherInformation = 1./Variance.^2;

% Compute embedding change probabilities and execute embedding
FI = FisherInformation(:)';
    
% Ternary embedding change probabilities
beta = TernaryProbs(FI,Payload);

% Simulate embedding
Stego = Cover;
beta = 2 * beta;
r = rand(1,numel(Cover));
ModifPM1 = (r < beta);                % Cover elements to be modified by +-1
r = rand(1,numel(Cover));
Stego(ModifPM1) = Cover(ModifPM1) + 2*(round(r(ModifPM1))) - 1; % Modifying X by +-1
Stego(Stego>255) = 253;                    % Taking care of boundary cases
Stego(Stego<0)   = 2;
ChangeRate = sum(ModifPM1(:))/numel(Cover); % Computing the change rate
pChange = reshape(beta/2,size(Cover));

end

% Beginning of the supporting functions

% Estimation of the pixels' variance 
function EstimatedVariance = VarianceEstimation(Image)

Kernel = ones(3,3);
% Local sums
x1 = conv2(Image, Kernel, 'same');
% Local quadratic sums
x2 = conv2(Image.^2, Kernel, 'same');
% Number of matrix elements in each square region
R = conv2(ones(size(Image)), Kernel, 'same');
% Local variance
EstimatedVariance = x2./R-(x1./R).^2;

end

% Computing the embedding change probabilities
function [beta] = TernaryProbs(FI,alpha)

load('ixlnx3.mat');

% Absolute payload in nats
payload = alpha * length(FI) * log(2);

% Initial search interval for lambda
[L, R] = deal (10^3, 10^6);

fL = h_tern(1./invxlnx3_fast(L*FI,ixlnx3)) - payload;
fR = h_tern(1./invxlnx3_fast(R*FI,ixlnx3)) - payload;
% If the range [L,R] does not cover alpha enlarge the search interval
while fL*fR > 0
    if fL > 0
        R = 2*R;
        fR = h_tern(1./invxlnx3_fast(R*FI,ixlnx3)) - payload;
    else
        L = L/2;
        fL = h_tern(1./invxlnx3_fast(L*FI,ixlnx3)) - payload;
    end
end

% Search for the labmda in the specified interval
[i, fM, TM] = deal(0, 1, zeros(60,2));
while (abs(fM)>0.0001 && i<60)
    M = (L+R)/2;
    fM = h_tern(1./invxlnx3_fast(M*FI,ixlnx3)) - payload;
    if fL*fM < 0, R = M; fR = fM;
    else          L = M; fL = fM; end
    i = i + 1;
    TM(i,:) = [fM,M];
end
if (i==60)
    M = TM(find(abs(TM(:,1)) == min(abs(TM(:,1))),1,'first'),2);
end
% Compute beta using the found lambda
beta = 1./invxlnx3_fast(M*FI,ixlnx3);

end

% Fast solver of y = x*log(x-2) paralellized over all pixels
function x = invxlnx3_fast(y,f)

i_large = y>1000;
i_small = y<=1000;

iyL = floor(y(i_small)/0.01)+1;
iyR = iyL + 1;
iyR(iyR>100001) = 100001;

x = zeros(size(y));
x(i_small) = f(iyL) + (y(i_small)-(iyL-1)*0.01).*(f(iyR)-f(iyL));

z = y(i_large)./log(y(i_large)-2);
for j = 1 : 20
    z = y(i_large)./log(z-2);
end
x(i_large) = z;

end

% Ternary entropy function expressed in nats
function Ht = h_tern(Probs)

p0 = 1-2*Probs;
P = [p0(:);Probs(:);Probs(:)];
H = -(P .* log(P));
H((P<eps)) = 0;
Ht = nansum(H);

end
