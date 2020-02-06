function [Stego, pChange, ChangeRate] = MiPODEx (cover, payload)
cover = single( imread(cover) );

% Compute Variance and do the flooring for numerical stability
WienerResidual = cover - wiener2(cover,[2,2]);
Variance = VarianceEstimationDCT2D(WienerResidual,9,9);
Variance(Variance< 0.01) = 0.01;

% Compute Fisher information and smooth it
% FisherInformation = 1./Variance.^2;
% FisherInformation = imfilter(FisherInformation,fspecial('average',7),'symmetric');
%% ÐÞ¸ÄÎª
FisherInformation = 1./Variance.^3;
FisherInformation = ordfilt2(FisherInformation, 7*7, true(7), 'symmetric');


% Compute embedding change probabilities and execute embedding
FI = FisherInformation(:)';
    
% Ternary embedding change probabilities
beta = TernaryProbs(FI,payload);

% Simulate embedding
Stego = cover;
beta = 2 * beta;
r = rand(1,numel(cover));
ModifPM1 = (r < beta);                % Cover elements to be modified by +-1
r = rand(1,numel(cover));
Stego(ModifPM1) = cover(ModifPM1) + 2*(round(r(ModifPM1))) - 1; % Modifying X by +-1
Stego(Stego>255) = 253;                    % Taking care of boundary cases
Stego(Stego<0)   = 2;
ChangeRate = sum(ModifPM1(:))/numel(cover); % Computing the change rate
pChange = reshape(beta/2,size(cover));

end

% Beginning of the supporting functions

% Estimation of the pixels' variance based on a 2D-DCT (trigonometric polynomial) model
function EstimatedVariance = VarianceEstimationDCT2D(Image, BlockSize, Degree)
% verifying the integrity of input arguments
if ~mod(BlockSize,2)
    error('The block dimensions should be odd!!');
end
if (Degree > BlockSize)
    error('Number of basis vectors exceeds block dimension!!');
end

% number of parameters per block
q = Degree*(Degree+1)/2;

% Build G matirx
BaseMat = zeros(BlockSize);BaseMat(1,1) = 1;
G = zeros(BlockSize^2,q);
k = 1;
for xShift = 1 : Degree
    for yShift = 1 : (Degree - xShift + 1)
        G(:,k) = reshape(idct2(circshift(BaseMat,[xShift-1 yShift-1])),BlockSize^2,1);
        k=k+1;
    end
end

% Estimate the variance
PadSize = floor(BlockSize/2*[1 1]);
I2C = im2col(padarray(Image,PadSize,'symmetric'),BlockSize*[1 1]);
PGorth = eye(BlockSize^2) - (G*((G'*G)\G'));
EstimatedVariance = reshape(sum(( PGorth * I2C ).^2)/(BlockSize^2 - q),size(Image));
end

% Computing the embedding change probabilities
function [beta] = TernaryProbs(FI,alpha)

load('D:\Matlab2017b\bin\embedAlg\MiPOD\ixlnx3.mat');

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
