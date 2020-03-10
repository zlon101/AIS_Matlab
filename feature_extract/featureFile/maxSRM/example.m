%               EXAMPLE - USING MEX "Spatial Rich Model"
%
% -------------------------------------------------------------------------
% Copyright (c) 2012 DDE Lab, Binghamton University, NY.
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
%          May 2012, http://dde.binghamton.edu
% -------------------------------------------------------------------------
% [1] Rich Models for Steganalysis of Digital Images, J. Fridrich and J.
% Kodovsky, IEEE Transactions on Information Forensics and Security, 2011.
% Under review. 
% http://dde.binghamton.edu/kodovsky/pdf/TIFS2012-SRM.pdf
% -------------------------------------------------------------------------
clc; clear all;

% Specify all images for extraction
I = imread('E:\astego\Images\BOSS_ALL\1.pgm');
MAP = ones(size(I));

%% --------------------
% SRM extraction by MEX 
% ---------------------
fprintf('maxSRM extraction');
MEXstart = tic;

%% Run default SRM extraction
F = maxSRM(I, MAP);

%% SRM MEX extraction can be configured in following way:
% (Using also default values - for detailed information see [1])
% Any or all the following settings might be included

%{
config.T = int32(2);
config.order = int32(4);
config.merge_spams = logical(true);
config.symm_sign = logical(true);
config.symm_reverse = logical(true);
config.symm_minmax = logical(true);
config.eraseLSB = logical(false);
config.parity = logical(false);

% F = SRM(ImageSet, config);
%}
        
MEXend = toc(MEXstart);
fprintf(' - DONE');
fprintf('\n\nThe image was extracted in %.2f seconds.\n', MEXend);
Ss = fieldnames(F);
fprintf('\n"F" contains %d submodel feature matrices (number of images x feature dimension): \n', numel(Ss));
for Sid = 1:length(Ss)
    Fsingle = eval(['F.' Ss{Sid}]);
    fprintf('   F.%s : %d x %d\n', Ss{Sid}, size(Fsingle, 1), size(Fsingle, 2));
end