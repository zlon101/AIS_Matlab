clear; clc;

addpath('costs')
addpath('distortion')
addpath('embedding')

payload = 0.4;

% load cover
cover = imread( '1013.pgm' );
% prepare costs
cover_rho = HILL( cover );
% if you want, this is how you change the costs from base distortion
% function costs into additive approximation of new costs IF you are using
% the weight parameter alpha = 0.4
cost = imfilter( cover_rho, 0.4/2*[0 1 0; 1 4 1; 0 1 0], 'conv', ...
    'symmetric' );

% if you want to start from cover
stego_start = cover;
% alternatively you can start at the additive approximation
rhoM1 = cost; rhoM1( cover ==   0 ) = 10^10;
rho0  = zeros( size( cost ) );
rhoP1 = cost; rhoP1( cover == 255 ) = 10^10;
stego_start = uint8( EmbeddingSimulator( cover, rhoM1, rho0, rhoP1, ...
    payload ) );
% the algorithm still needs the cover for the distortion function, the role
% of the stego_start is to allow to continue computing sweeps when you
% already computed several and decided it is not enough

% clean up so the variable names after rhe next steps are not confusing
clear cost rhoM1 rho0 rhoP1

% how many sweeps do you want the Gibbs construction to run
sweeps = 5;

%% the construction itself
% IF YOU CAN, compile a mex file, it runs much faster
[stego_all, distortion_all, costs_all, lambda_all, grid] = ...
    Gibbs( cover, cover_rho, stego_start, sweeps, payload );

%% stego_all contains all stego images, first one is the stego_start
% for the purposes of imwrite all images are in uint8
cover = double( cover );
stego_start = double( stego_start );
stego_all = double( stego_all );
figure
subplot(1,3,1)
imshow( cover, [0 255] );
title('cover')
subplot(1,3,2)
imshow( stego_all(:,:,1) - cover, [] );
title('additive approximation')
subplot(1,3,3)
imshow( stego_all(:,:,2) - cover, [] );
title('1st sweep')

%% distortion_all holds all distortions, the first one being the distortion
% of the stego_start with respect to cover
figure
plot( 0:sweeps, distortion_all )
xlabel('Sweeps')
ylabel('Distortion')

%% costs_all.cover_rho has a copy of the original costs
% costs_all.rhoM1, .rho0 and .rhoP1 have the hold the costs of each pixel
% after each sweep, for all three the costs_all.rhoXX(:,:,1) are zeros and
% are stored only so the indexing in sweeps is consistent through all saved
% data
figure
rho = zeros( size( cover_rho ) );
rhoM1 = costs_all.rhoM1(:,:,2);
imshow( nowet( costs_all.rhoM1(:,:,2) ), [] );
% nowet is a function for imshow, changes costs larger to 10^7 to the
% largest known cost smaller than 10^7 so the image of costs is not 
% oversaturated with wet costs
title('Costs of embedding -1 used in the 1st sweep')

%% lambda_all hold all lambdas used in the embedding simulator, rows are
% sweeps, columns are lattices, once again, first row are 0s to keep
% indexing consistent
figure
hold on
plot(1:sweeps, lambda_all(2:end,1), 'b')
plot(1:sweeps, lambda_all(2:end,2), 'r')
hold off
xlabel('Sweep')
legend('\lambda for the 1st lattice', '\lambda for the 2nd lattice')

%% last but not least, grid tells you which pixel belongs to which lattice
% and also the visiting order -- the lattice 1 is always visited first
figure
imshow( grid == 1, [] )
title('Pixels in the 1st lattice')