function stego= embedAlgCZL(cover,payload,Am)
% 失真函数设计
%%
% sharpImg= sharpen(single(imread(cover)), str2double(Am));
% [rhoP1,rhoM1]= CostHUGO_mex( sharpImg );
% -------------------------------------------------------------------------
cover= single(imread(cover));
[rhoP1,rhoM1]= CMeanFilterT13( sharpen(cover,Am) );

stego= EmbeddingSimulator(cover, single(rhoP1), single(rhoM1), payload*numel(cover), false);
end

%% HUGO
function [rhoP1,rhoM1] = Cost_CZL(cover)
% HUGO 代价函数
% 返回+1 -1 的代价
%% 
cover = single(cover);
wetCost = 10^8;

% create mirror padded cover image
padSize = double(3);
cPadded = padarray(cover, [padSize,padSize], 'symmetric');
% create residuals
rezH = (cPadded(:, 1:end-1) - cPadded(:, 2:end)).^2;
rezV = (cPadded(1:end-1, :) - cPadded(2:end, :)).^2;
% rezD = cPadded(1:end-1, 1:end-1) - cPadded(2:end, 2:end);
% rezMD = cPadded(1:end-1, 2:end) - cPadded(2:end, 1:end-1);

rhoM1 = zeros(size(cover),'single');  % declare cost of -1 change           
rhoP1 = zeros(size(cover),'single');  % declare cost of +1 change
T =3; a=2; % 滤波器阶数 & 权重, T=3,5,7
T=(T-1)*0.5;
cH=1; cV=1;
for row=1:size(cover, 1)
  r=row+3;
  for col=1:size(cover, 2)
    c=col+3; % padSize=3;
    % T=3
    arrH= [rezH(r-1, c-T:c+T-1);
           rezH(r,   c-T:c+T-1).*a;
           rezH(r+1, c-T:c+T-1)];
    arrV= [rezV(r-T:r+T-1, c-1);
           rezV(r-T:r+T-1, c).*a;
           rezV(r-T:r+T-1, c+1)];     
    % CZL_8
    %DCover = norm( reshape([arrDH;a*rezH],[],1) )+ norm( reshape([arrDV;a*rezV],[],1) );
    %DCover=0.5*DCover;
    %rhoP1(row,col)=1./(DCover+1e-20);
    
    % CZL7
    %DCover = norm( reshape([arrDH;a*rezH],[],1) );
    %rhoP1(row,col)= 1./(DCover+1e-20);
  end
end

% 平滑滤波
L= ones(9);
rhoP1= imfilter(rhoP1, L,'symmetric','conv','same')./sum(L(:));
% rhoP1 = ordfilt2(rhoP1,81,true(9),'symmetric');

rhoM1= rhoP1;

rhoM1(rhoM1>wetCost) = wetCost;
rhoP1(rhoP1>wetCost) = wetCost;
rhoP1(cover == 255) = wetCost;
rhoM1(cover == 0) = wetCost;

%{
t =1; a=2; % 滤波器阶数 & 权重
for row=1:size(cover, 1)
  for col=1:size(cover, 2)
    % Horizontal
    arrDH=[CrezH(row+2, col:col+t);
           CrezH(row+4, col:col+t)];
    rezH = CrezH(row+3, col:col+t);
%}
end

%{
%% SUNWD
function [rhoP1,rhoM1] = CostUNWDOpt(coverImg)
if(ischar(coverImg))
	coverImg = imread(coverImg); 
end
cover = single(coverImg);
[k,l] = size(cover);
sgm = 1;
wetCost = 10^8;
%% Get 2D wavelet filters - Daubechies 8
% 1D high pass decomposition filter
hpdf = [-0.0544158422, 0.3128715909, -0.6756307363, 0.5853546837, 0.0158291053, -0.2840155430, -0.0004724846, 0.1287474266, 0.0173693010, -0.0440882539, ...
        -0.0139810279, 0.0087460940, 0.0048703530, -0.0003917404, -0.0006754494, -0.0001174768];
% 1D low pass decomposition filter
lpdf = (-1).^(0:numel(hpdf)-1).*fliplr(hpdf);
% construction of 2D wavelet filters
F{1} = lpdf'*hpdf;
F{2} = hpdf'*lpdf;
F{3} = hpdf'*hpdf;

% add padding
padSize = max([size(F{1})'; size(F{2})'; size(F{3})']);
coverPadded = padarray(cover, [padSize padSize], 'symmetric');

xi = cell(3, 1);
for fIndex = 1:3
  % compute residual
  R = conv2(coverPadded, F{fIndex}, 'same');
  % compute suitability
  xi{fIndex} = conv2(1./(abs(R)+sgm), rot90(abs(F{fIndex}),2), 'same');
  % correct the suitability shift if filter size is even
  if mod(size(F{fIndex}, 1), 2) == 0
    xi{fIndex} = circshift(xi{fIndex}, [1, 0]);
  end
  if mod(size(F{fIndex}, 2), 2) == 0
    xi{fIndex} = circshift(xi{fIndex}, [0, 1]); 
  end
  % remove padding
  xi{fIndex} = xi{fIndex}(((size(xi{fIndex}, 1)-k)/2)+1:end-((size(xi{fIndex}, 1)-k)/2), ((size(xi{fIndex}, 2)-l)/2)+1:end-((size(xi{fIndex}, 2)-l)/2));
end

% compute embedding costs \rho
rho = xi{1} + xi{2} + xi{3};
L = ones(5);
rho = imfilter(rho, L,'symmetric','conv','same')./sum(L(:));

% adjust embedding costs
rho(rho > wetCost) = wetCost; % threshold on the costs
rho(isnan(rho)) = wetCost; % if all xi{} are zero threshold the cost
rhoP1 = rho;  % +1 的代价
rhoM1 = rho;
rhoP1(cover==255) = wetCost; % do not embed +1 if the pixel has max value
rhoM1(cover==0) = wetCost; % do not embed -1 if the pixel has min value
end

%% WOW
function [rhoP1,rhoM1] = CostWOWOpt(cover)
%% Get 2D wavelet filters - Daubechies 8
% 1D high pass decomposition filter
hpdf = [-0.0544158422, 0.3128715909, -0.6756307363, 0.5853546837, 0.0158291053, -0.2840155430, -0.0004724846, 0.1287474266, 0.0173693010, -0.0440882539, ...
        -0.0139810279, 0.0087460940, 0.0048703530, -0.0003917404, -0.0006754494, -0.0001174768];
% 1D low pass decomposition filter
lpdf = (-1).^(0:numel(hpdf)-1).*fliplr(hpdf);
% construction of 2D wavelet filters
F{1} = lpdf'*hpdf;
F{2} = hpdf'*lpdf;
F{3} = hpdf'*hpdf;

%% Get embedding costs
% inicialization
cover = double(cover);
params.p = -1;  % holder norm parameter
p = params.p;
wetCost = 10^10;
sizeCover = size(cover);

% add padding
padSize = max([size(F{1})'; size(F{2})'; size(F{3})']);
coverPadded = padarray(cover, [padSize padSize], 'symmetric');

% compute directional residual and suitability \xi for each filter
xi = cell(3, 1);
for fIndex = 1:3
  % compute residual
  R = conv2(coverPadded, F{fIndex}, 'same');

  % compute suitability
  xi{fIndex} = conv2(abs(R), rot90(abs(F{fIndex}), 2), 'same');
  % correct the suitability shift if filter size is even
  if mod(size(F{fIndex},1), 2) == 0
      xi{fIndex} = circshift(xi{fIndex}, [1, 0]); 
  end
  if mod(size(F{fIndex}, 2), 2) == 0
      xi{fIndex} = circshift(xi{fIndex}, [0, 1]); 
  end
  % remove padding
  xi{fIndex} = xi{fIndex}(((size(xi{fIndex}, 1)-sizeCover(1))/2)+1:end-((size(xi{fIndex}, 1)-sizeCover(1))/2), ((size(xi{fIndex}, 2)-sizeCover(2))/2)+1:end-((size(xi{fIndex}, 2)-sizeCover(2))/2));
end

% compute embedding costs \rho
rho = ( (xi{1}.^p) + (xi{2}.^p) + (xi{3}.^p) ) .^ (-1/p);
%% 平滑
L2 = ones(5);
rho = imfilter(rho, L2,'symmetric', 'conv', 'same')./sum(L2(:));
% adjust embedding costs
rho(rho > wetCost) = wetCost; % threshold on the costs
rho(isnan(rho)) = wetCost; % if all xi{} are zero threshold the cost
rhoP1 = rho;
rhoM1 = rho;
rhoP1(cover==255) = wetCost; % do not embed +1 if the pixel has max value
rhoM1(cover==0) = wetCost; % do not embed -1 if the pixel has min value
end

%% HILL
function [rhoP1,rhoM1]=HILLOpt(I)
HFilter = [-1, 2, -2, 2,-1;
           2,-6,  8,-6, 2;
          -2, 8,-12, 8,-2;
           2,-6,  8,-6, 2;
          -1, 2, -2, 2,-1];
HF2 = [1,2,1;2,-12,2;1,2,1];
KB = [-1,2,-1; 2,-4,2; -1,2,-1];
L1= ones(3);
L2 = ones(15);

Y1 = imfilter(I, KB,'symmetric','same');
Y2 = imfilter(abs(Y1), L1,'symmetric','same')./sum(L1(:));
% 滤波平滑
Cost = imfilter((Y2+1e-20).^-1, L2,'symmetric','same')./sum(L2(:));
Cost = ordfilt2(Cost, 121, true(11), 'symmetric');

wetCost = 10^8;
Cost(Cost > wetCost) = wetCost;
Cost(isnan(Cost)) = wetCost;
rhoP1 = Cost;  % +1 的代价
rhoM1 = Cost;
rhoP1(I==255) = wetCost;
rhoM1(I==0) = wetCost;
end
%}