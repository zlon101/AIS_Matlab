function F = CSR( image_path )
% -------------------------------------------------------------------------
% Copyright (c) 2014 DDE Lab, Binghamton University, NY.
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
% Authors: Tomas Denemark, Jessica Fridrich, Vojtech Holub
% -------------------------------------------------------------------------
% Contact: tdenema1@binghamton.edu
%          fridrich@binghamton.edu
%          vojtech_holub@yahoo.com
%          January 2014, http://dde.binghamton.edu
% -------------------------------------------------------------------------
% Extracts CSR features described in [1] from a single image using default
% settings using all six filters. Total dimensionality of the features
% vector is 1163. This script uses a slightly moddified S-UNIWARD [2][3]
% embedding simulator that is attached.
% -------------------------------------------------------------------------
% F = CSR( image_path ) extracts the features from image at 'image_path'
% with default settings
%--------------------------------------------------------------------------
% Input format:
%   image_path     - string with path to the input image
%
% Output format:
%   F              - 1183 features extracted from the image
% -------------------------------------------------------------------------
% [1] Further Study on the Security of S-UNIWARD, T. Denemark and
%     J. Fridrich and V. Holub , Proc. SPIE, Electronic Imaging, Media
%     Watermarking, Security, and Forensics, 2014
%     http://dde.binghamton.edu/tomasD/pdf/ ...
%     SPIE14_Further_Study_on_Security_of_S-UNIWARD.pdf
% [2] Universal Distortion Function for Steganography in an Arbitrary 
%     Domain V. Holub, J. Fridrich and T. Denemark, EURASIP Journal on 
%     Information Security 2014(1) (Section: SI: Revised Selected Papers 
%     of ACM IH and MMS 2013) 
% [3] Digital Image Steganography Using Universal Distortion, V. Holub and
%     J. Fridrich, Proc. ACM Workshop on Information Hiding and Multimedia
%     Security, 2013
%     http://dde.binghamton.edu/vholub/pdf/ ...
%     ACMIH2013_Digital_Image_Steganography_Using_Universal_Distortion.pdf
% -------------------------------------------------------------------------

payload_bar = 0.4;      % testing payload
sigma_bar = 10*eps;     % testing stabilization constant
Ts = 0.05;              % threshold for 'small' probability class
TL = 0.06;              % threshold for 'Large' probability class
filters = {'1st1D', '1st2D', '2nd1D', '2nd2D', '3rd1D'};
                        % cell aray of filters to be used

FCell = cell(1, 1);     % dummy variable

X = double( imread( image_path ));
                        % load the image
[~, P] = S_UNIWARD(X, payload_bar, sigma_bar);
                        % extract the embedding probabilities
L = P > TL; s = P < Ts; % 'small' and 'Large' classes

% Extract CSR features
for fi = 1:length(filters)
    filter = filters{fi};
    switch filter
        case '1st1D'
            f = @f1st1D;
            T = 10;
        case '1st2D'
            f = @f1st2D;
            T = 3;
        case '1st3D'
            f = @f1st3D;
            T = 2;
        case '2nd1D'
            f = @f2nd1D;
            T = 10;
        case '2nd2D'
            f = @f2nd2D;
            T = 3;
        case '3rd1D'
            f = @f3rd1D;
            T = 10;
    end
    FCell{fi} = f(X, L, s, T);
                        % extraction of features for each filter
end

F = [];                 % variable declaration
for i = 1:length(filters)
    F = [F FCell{i}];   % result collection
end

end

% -------------------------------------------------------------------------
% Individiual feature extractors
% -------------------------------------------------------------------------

function f = f1st1D(X, L, s, T)

X = [X, X'];            % Adding vertical direction to the image
L = [L, L'];
s = [s, s'];

Centr =  X(:,1:end-1) - X(:,2:end);
                        % 1st-order CENTRAL horizontal residual

Centr =  min( max(Centr, -T), T);
                        % Truncate

% PURE classes
% [s s] class
inc = incidence('ss', s, L);
                        % Incidence array
hss = reshape(cooc(T, Centr(inc)), 1,[]);
                        % Features as a 1D co-occurence

% [L L] class
inc = incidence('LL', s, L);
                        % Incidence array
hLL = reshape(cooc(T, Centr(inc)), 1,[]);
                        % Features as a 1D co-occurence

% MIXED classes
% [s L] and [L ] classes
inc = incidence('sL', s, L); h1 = cooc(T, Centr(inc));
                        % Features for the first class
inc = incidence('Ls', s, L); h2 = cooc(T, Centr(inc));
                        % Features for the second class
hsL = combf(h1, h2); % Combine the features

% Forming of the feature vector
f = [hss hLL hsL ];
end
function f = f1st2D(X, L, s, T)

X = [X, X'];            % Adding vertical direction to the image
L = [L, L'];
s = [s, s'];

Left =  X(:,1:end-2) - X(:,2:end-1);
                        % 1st-order LEFT  horizontal residual
Right = X(:,2:end-1) - X(:,3:end);
                        % 1st-order RIGHT horizontal residual

Left =  min( max(Left,  -T), T);
                        % Truncate
Right = min( max(Right, -T), T);
                        % Truncate

% PURE classes
% [s s s] class
inc = incidence('sss', s, L);
                        % Incidence array
hsss = reshape(cooc(T, [Left(inc), Right(inc)]), 1,[]);
                        % Features as a 2D co-occurence

% [s L s] class
inc = incidence('sLs', s, L);
                        % Incidence array
hsLs = reshape(cooc(T, [Left(inc), Right(inc)]), 1,[]);
                        % Features as a 2D co-occurence

% [L s L] class
inc = incidence('LsL', s, L);
                        % Incidence array
hLsL = reshape(cooc(T, [Left(inc), Right(inc)]), 1,[]);
                        % Features as a 2D co-occurence

% [L L L] class
inc = incidence('LLL', s, L);
                        % Incidence array
hLLL = reshape(cooc(T, [Left(inc), Right(inc)]), 1,[]);
                        % Features as a 2D co-occurence

% MIXED classes
% [s s L] and [L s s] classes
inc = incidence('ssL', s, L); h1 = cooc(T, [Left(inc), Right(inc)]);
                        % Features for the first class
inc = incidence('Lss', s, L); h2 = cooc(T, [Left(inc), Right(inc)]);
                        % Features for the second class
hssL = combf(h1, h2);
                        % Combine the features

% [s L L] and [L L s] classes
inc = incidence('sLL', s, L); h1 = cooc(T, [Left(inc), Right(inc)]);
inc = incidence('LLs', s, L); h2 = cooc(T, [Left(inc), Right(inc)]);
hsLL = combf(h1, h2);

% Forming of the feature vector
f = [hsss hsLs hLsL hLLL hssL hsLL ];

end
function f = f2nd1D(X, L, s, T)

X = [X, X'];            % Adding vertical direction to the image
L = [L, L'];
s = [s, s'];

Centr =  X(:,1:end-2) - 2*X(:,2:end-1) + X(:,3:end);
                        % 1st-order LEFT  horizontal residual
Centr = min( max(Centr, -T), T);
                        % Truncate

% PURE classes
% [s s s] class
inc = incidence('sss', s, L);
                        % Incidence array
hsss = reshape(cooc(T, Centr(inc)), 1,[]);
                        % Features as a 2D co-occurence

% [s L s] class
inc = incidence('sLs', s, L);
                        % Incidence array
hsLs = reshape(cooc(T, Centr(inc)), 1,[]);
                        % Features as a 2D co-occurence

% [L s L] class
inc = incidence('LsL', s, L);
                        % Incidence array
hLsL = reshape(cooc(T, Centr(inc)), 1,[]);
                        % Features as a 2D co-occurence

% [L L L] class
inc = incidence('LLL', s, L);
                        % Incidence array
hLLL = reshape(cooc(T, Centr(inc)), 1,[]);
                        % Features as a 2D co-occurence

% MIXED classes
% [s s L] and [L s s] classes
inc = incidence('ssL', s, L); h1 = cooc(T, Centr(inc));
                        % Features for the first class
inc = incidence('Lss', s, L); h2 = cooc(T, Centr(inc));
                        % Features for the second class
hssL = combf(h1, h2);
                        % Combine the features

% [s L L] and [L L s] classes
inc = incidence('sLL', s, L); h1 = cooc(T, Centr(inc));
inc = incidence('LLs', s, L); h2 = cooc(T, Centr(inc));
hsLL = combf(h1, h2);

% Forming of the feature vector
f = [hsss hsLs hLsL hLLL hssL hsLL ];
        
end
function f = f2nd2D(X, L, s, T)

X = [X,X'];             % Adding vertical direction to the image
L = [L,L'];
s = [s,s'];

Left =  X(:,1:end-3) - 2*X(:,2:end-2) + X(:,3:end-1);
                        % 2nd-order LEFT  horizontal residual
Right = X(:,2:end-2) - 2*X(:,3:end-1) + X(:,4:end);
                        % 2nd-order RIGHT horizontal residual

Left =  min( max(Left,  -T), T);
                        % Truncate
Right = min( max(Right, -T), T);
                        % Truncate

% PURE classes
% [s s s s] class
inc = incidence('ssss', s, L);
                        % Incidence array
hssss = reshape(cooc(T, [Left(inc), Right(inc)]), 1,[]);
                        % Features as a 2D co-occurence

% [s L L s] class
inc = incidence('sLLs', s, L);
                        % Incidence array
hsLLs = reshape(cooc(T, [Left(inc), Right(inc)]), 1,[]);
                        % Features as a 2D co-occurence

% [L s s L] class
inc = incidence('LssL', s, L);
                        % Incidence array
hLssL = reshape(cooc(T, [Left(inc), Right(inc)]), 1,[]);
                        % Features as a 2D co-occurence

% [L L L L] class
inc = incidence('LLLL', s, L);
                        % Incidence array
hLLLL = reshape(cooc(T, [Left(inc), Right(inc)]), 1,[]);
                        % Features as a 2D co-occurence

% MIXED classes
% [s s s L] and [L s s s] classes
inc = incidence('sssL', s, L); h1 = cooc(T, [Left(inc), Right(inc)]);
                        % Features for the first class
inc = incidence('Lsss', s, L); h2 = cooc(T, [Left(inc), Right(inc)]);
                        % Features for the second class
hsssL = combf(h1, h2);
                        % Combine the features

% [s s L s] and [s L s s] classes
inc = incidence('ssLs', s, L); h1 = cooc(T, [Left(inc), Right(inc)]);
inc = incidence('sLss', s, L); h2 = cooc(T, [Left(inc), Right(inc)]);
hssLs = combf(h1, h2);

% [s s L L] and [L L s s] classes
inc = incidence('ssLL', s, L); h1 = cooc(T, [Left(inc), Right(inc)]);
inc = incidence('LLss', s, L); h2 = cooc(T, [Left(inc), Right(inc)]);
hssLL = combf(h1, h2);

% [s L s L] and [L s L s] classes
inc = incidence('sLsL', s, L); h1 = cooc(T, [Left(inc), Right(inc)]);
inc = incidence('LsLs', s, L); h2 = cooc(T, [Left(inc), Right(inc)]);
hsLsL = combf(h1, h2);

% [s L L L] and [L L L s] classes
inc = incidence('sLLL', s, L); h1 = cooc(T, [Left(inc), Right(inc)]);
inc = incidence('LLLs', s, L); h2 = cooc(T, [Left(inc), Right(inc)]);
hsLLL = combf(h1, h2);

% [L s L L] and [L L s L] classes
inc = incidence('LsLL', s, L); h1 = cooc(T, [Left(inc), Right(inc)]);
inc = incidence('LLsL', s, L); h2 = cooc(T, [Left(inc), Right(inc)]);
hLsLL = combf(h1, h2);

% Forming of the feature vector
f = [hssss hsLLs hLssL hLLLL hsssL hssLs hssLL hsLsL hsLLL hLsLL ];

end
function f = f3rd1D(X, L, s, T)

X = [X, X'];            % Adding vertical direction to the image
L = [L, L'];
s = [s, s'];

Centr =  -X(:,1:end-3) + 3*X(:,2:end-2) - 3*X(:,3:end-1) + X(:,4:end);
                        % 3rd-order LEFT  horizontal residual
Centr =  min( max(Centr,  -T), T);
                        % Truncate

% PURE classes
% [s s s s] class
inc = incidence('ssss', s, L);
                        % Incidence array
hssss = reshape(cooc(T, Centr(inc)), 1,[]);
                        % Features as a 2D co-occurence

% [s L L s] class
inc = incidence('sLLs', s, L);
                        % Incidence array
hsLLs = reshape(cooc(T, Centr(inc)), 1,[]);
                        % Features as a 2D co-occurence

% [L s s L] class
inc = incidence('LssL', s, L);
                        % Incidence array
hLssL = reshape(cooc(T, Centr(inc)), 1,[]);
                        % Features as a 2D co-occurence

% [L L L L] class
inc = incidence('LLLL', s, L);
                        % Incidence array
hLLLL = reshape(cooc(T, Centr(inc)), 1,[]);
                        % Features as a 2D co-occurence

% MIXED classes
% [s s s L] and [L s s s] classes
inc = incidence('sssL', s, L); h1 = cooc(T, Centr(inc));
                        % Features for the first class
inc = incidence('Lsss', s, L); h2 = cooc(T, Centr(inc));
                        % Features for the second class
hsssL = combf(h1, h2);  % Combine the features

% [s s L s] and [s L s s] classes
inc = incidence('ssLs', s, L); h1 = cooc(T, Centr(inc));
inc = incidence('sLss', s, L); h2 = cooc(T, Centr(inc));
hssLs = combf(h1, h2);

% [s s L L] and [L L s s] classes
inc = incidence('ssLL', s, L); h1 = cooc(T, Centr(inc));
inc = incidence('LLss', s, L); h2 = cooc(T, Centr(inc));
hssLL = combf(h1, h2);

% [s L s L] and [L s L s] classes
inc = incidence('sLsL', s, L); h1 = cooc(T, Centr(inc));
inc = incidence('LsLs', s, L); h2 = cooc(T, Centr(inc));
hsLsL = combf(h1, h2);

% [s L L L] and [L L L s] classes
inc = incidence('sLLL', s, L); h1 = cooc(T, Centr(inc));
inc = incidence('LLLs', s, L); h2 = cooc(T, Centr(inc));
hsLLL = combf(h1, h2);

% [L s L L] and [L L s L] classes
inc = incidence('LsLL', s, L); h1 = cooc(T, Centr(inc));
inc = incidence('LLsL', s, L); h2 = cooc(T, Centr(inc));
hLsLL = combf(h1, h2);

% Forming of the feature vector
f = [hssss hsLLs hLssL hLLLL hsssL hssLs hssLL hsLsL hsLLL hLsLL ];

end

% -------------------------------------------------------------------------
% Auxiliary functions
% -------------------------------------------------------------------------

function inc = incidence(class, s, L)

n = length(class);
inc = zeros(size(s,1), size(s,2)-n+1);
                        % declaration of inc

for i = 1:n             % for every letter in the class definition
    switch class(i)
        case 's'        % 'small' probabilities
            inc = inc + s(:, i:end-n+i);
        case 'L'        % 'large' probabilities
            inc = inc + L(:, i:end-n+i);
        otherwise
            error();
    end
end

inc = inc == n;         % convert to logical, find only the correct classes

end
function h = combf(h1, h2)

[D(1), D(2), D(3)] = size(h1);
D = sum(D > 1);        % the dimension of the feature set
switch D
    case 1
        h = h1 + h2;
    case 2
        h = h1 + h2';
    case 3
        h = zeros(size(h1));
        for i = 1:size(h,1)
            for j = 1:size(h,2)
                for k = 1:size(h,2)
                    h(i,j,k) = h1(i,j,k) + h2(i,j,k);
                end
            end
        end
end
h = h(:)';

end

function f = cooc(T, v)
% f = cooc(T, v)
% performs nD co-occurence
% INPUTS:
% T -- threshold
% v -- vector
% OUTPUTS:
% f -- vector or field of featues

order = size(v,2); % order of the co-occurence

switch order
    case 1
        f = hist(v, -T:T);
    case 2
        L = v(:,1); R = v(:,2);
        f = zeros(2*T+1, 2*T+1);
        for i = -T : T
            R2 = R(L(:)==i);
            for j = -T : T
                f(i+T+1,j+T+1) = sum(R2(:)==j);
            end
        end
    case 3
        L = v(:,1); C = v(:,2); R = v(:,3);
        f = zeros(2*T+1, 2*T+1, 2*T+1);
        for i = -T : T
            C2 = C(L(:)==i);
            for j = -T : T
                R2 = R(C2(:)==j);
                for k = -T:T
                    f(i+T+1,j+T+1,k+T+1) = sum(R2(:)==k);
                end
            end
        end
end

end
