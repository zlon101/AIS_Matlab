function f = cchen972(IMAGE,QF)
% -------------------------------------------------------------------------
% Copyright (c) 2011 DDE Lab, Binghamton University, NY.
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
% Contact: jan@kodovsky.com | fridrich@binghamton.edu | November 2011
%          http://dde.binghamton.edu/download/feature_extractors
% -------------------------------------------------------------------------
% Extracts 486 inter- and intra-block DCT domain features proposed in [1].
% Additionally, this extractor enhaces these features by Cartesian
% Calibration as described in [3], using imwrite JPEG compressor. The
% dimensionality thus doubles to 2x486=972.
% -------------------------------------------------------------------------
% Input:  IMAGE ... path to the JPEG image
%         QF ...... JPEG quality factor of the reference image (we
%                   recommend to keep it the same as the original image QF)
% Output: f ....... extracted CC-CHEN feature vector
% -------------------------------------------------------------------------
% IMPORTANT: For accessing DCT coefficients of JPEG images, we use Phil
% Sallee's Matlab JPEG toolbox (function jpeg_read) available at
% http://philsallee.com/jpegtbx/
% -------------------------------------------------------------------------
% [1] C. Chen and Y. Q. Shi. "JPEG image steganalysis utilizing both
%     intrablock and interblock correlations," IEEE ISCAS, International
%     Symposium on Circuits and Systems, pages 3029â€?032, May 2008.
% [2] Y. Q. Shi, C. Chen, and W. Chen. "A Markov process based approach to
%     effective attacking JPEG steganography," In J. L. Camenisch, C. S.
%     Collberg, N. F. Johnson, and P. Sallee, editors, Information Hiding,
%     8th International Workshop, volume 4437 of Lecture Notes in Computer
%     Science, pages 249â€?64, Alexandria, VA, July 10â€?2, 2006.
%     Springer-Verlag, New York.
% [3] J. Kodovsky and J. Fridrich. "Calibration revisited." In J. Dittmann,
%     S. Craver, and J. Fridrich, editors, Proceedings of the 11th ACM
%     Multimedia & Security Workshop, pages 63â€?4, Princeton, NJ, September
%     7â€?, 2009.
% -------------------------------------------------------------------------

X = DCTPlane(IMAGE);
F1 = extractShi324(X); % intra-block part [2]
F2 = extract_interblock(X); % inter-block part [1]
f = [F1;F2];

% reference version of the features (for calibration)
X = DCTPlane_reference(IMAGE,QF);
F1ref = extractShi324(X); % intra-block reference part
F2ref = extract_interblock(X); % inter-block reference part
f = [f;F1ref;F2ref]; % append second half of features for Cartesian calibration [3]

function F = extract_interblock(A)
T = 4; % threshold
F = zeros(2*(2*T+1)^2,1); % inter-block part
for MODE = 2:64
    % inter-block part (for single mode)
    % obtain appropriate mode 2-d array (see [1])
    AA = PlaneToVecMode(A,MODE);
    % horizontal (2T+1)^2 features
    Ad = conv2(AA,[-1 1],'valid');   % <=> A(:,1:end-1) - A(:,2:end);
    Ad = max(min(Ad,T),-T);          % truncate to [-T,T]
    Fah = getTPM(Ad(:,1:end-1),Ad(:,2:end),T);         % add to the feature vector

    % vertical (2T+1)^2 features
    Ad = conv2(AA,[-1;1],'valid');    % <=> A(1:end-1,:) - A(2:end,:);
    Ad = max(min(Ad,T),-T);          % truncate to [-T,T]
    Fav = getTPM(Ad(1:end-1,:),Ad(2:end,:),T);       % add to the feature vector
    
    F = F + [Fah;Fav];                   % add to the final vector
end
F = F/63; % take average over all modes
function F = extractShi324(A)
% extract [2] features from DCTPlane A
A = abs(A);
T = 4; F = [];

%% horizontal (2T+1)^2 features
Ad = conv2(A,[-1 1],'valid');    % <=> A(:,1:end-1) - A(:,2:end);
Ad = max(min(Ad,T),-T);          % truncate to [-T,T]
A1 = Ad(:,1:end-1);
A2 = Ad(:,2:end);
F = [F;getTPM(A1,A2,T)];  % Markov chain

%% vertical (2T+1)^2 features
Ad = conv2(A,[-1;1],'valid');    % <=> A(1:end-1,:) - A(2:end,:);
Ad = max(min(Ad,T),-T);          % truncate to [-T,T]
A1 = Ad(1:end-1,:);
A2 = Ad(2:end,:);
F = [F;getTPM(A1,A2,T)];  % Markov chain

%% diagonal (2T+1)^2 features
Ad = conv2(A,[-1 0;0 1],'valid');% <=> A(1:end-1,1:end-1) - A(2:end,2:end);
Ad = max(min(Ad,T),-T);          % truncate to [-T,T]
A1 = Ad(1:end-1,1:end-1);
A2 = Ad(2:end,2:end);
F = [F;getTPM(A1,A2,T)]; % Markov chain

%% minor diagonal (2T+1)^2 features
Ad = conv2(A,[0 1;-1 0],'valid');% <=> A(2:end,1:end-1) - A(1:end-1,2:end);
Ad = max(min(Ad,T),-T);          % truncate to [-T,T]
A1 = Ad(2:end,1:end-1);
A2 = Ad(1:end-1,2:end);
F = [F;getTPM(A1,A2,T)]; % Markov chain
function F = getTPM(A1,A2,T)
% get transition probability matrix A1 --> A2, range -T..T
F = zeros(2*T+1);
dn = max(hist(A1(:),-T:T),1); % normalization factors
for i=-T:T
    FF = A2(A1==i); % filtered version
    for j=-T:T
        F(i+T+1,j+T+1) = nnz(FF==j)/dn(i+T+1);
    end
end
F = F(:);
function Mat=PlaneToVecMode(plane,MODE)
mask = reshape(1:64,8,8);
[i,j] = find(mask==MODE);
Mat = plane(i:8:end,j:8:end);
function Plane=DCTPlane(path)
% loads DCT Plane of the given JPEG image
jobj=jpeg_read(path); % Phil Sallee's MATLAB jpeg toolbox needed
Plane=jobj.coef_arrays{1};
function Plane=DCTPlane_reference(path,QF)
% obtain reference image DCT plane (see [3])
I = imread(path);      % decompress into spatial domain
I = I(5:end-4,5:end-4); % crop by 4x4 pixels
TMP = ['img_' num2str(round(rand()*1e7)) num2str(round(rand()*1e7)) '.jpg']; % temporary reference image
while exist(TMP,'file'), TMP = ['img_' num2str(round(rand()*1e7)) num2str(round(rand()*1e7)) '.jpg']; end
imwrite(I,TMP,'Quality',QF); % save as temporary jpeg image using imwrite
Plane = DCTPlane(TMP); % load DCT plane of the reference image
delete(TMP); % delete the temporary reference image