%               EXAMPLE - USING MEX "Projection Spatial Rich Model"
%
% -------------------------------------------------------------------------
% Copyright (c) 2013 DDE Lab, Binghamton University, NY.
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
% Author: Vojtech Holub
% -------------------------------------------------------------------------
% Contact: vojtech_holub@yahoo.com
%          jan@kodovsky.com
%          fridrich@binghamton.edu
%          November 2012, http://dde.binghamton.edu
% -------------------------------------------------------------------------
clc; close all;

% Specify all images for extraction
I = imread('1.pgm');

tstart = tic;
%% --------------------
% Run default PSRM extraction
% ---------------------
q = 1; % Use q=3 for detection of JPEGs in spatial domain

fprintf('PSRM extraction');
F = PSRM(I, q);

%% Results        
tend = toc(tstart);
fprintf(' - DONE');
fprintf('\n\nPSRM features extracted %.2f seconds\n', tend);
Ss = fieldnames(F);
fprintf('\n"F" contains %d submodel feature matrices (number of images x feature dimension): \n', numel(Ss));
for Sid = 1:length(Ss)
    Fsingle = eval(['F.' Ss{Sid}]);
    fprintf('   F.%s : %d x %d\n', Ss{Sid}, size(Fsingle, 1), size(Fsingle, 2));
end