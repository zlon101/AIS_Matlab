function F=mainDctrMex(ImageSet)
%               EXAMPLE - USING MEX "DCTR features"
%
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
% Author: Vojtech Holub
% -------------------------------------------------------------------------
% Contact: vojtech_holub@yahoo.com
%          fridrich@binghamton.edu
%          March 2014, http://dde.binghamton.edu
% -------------------------------------------------------------------------
%clc; clear all;

% Specify all images for extraction
%ImageSet = {fullfile('..', '..', 'image_dir', 'image.jpg')};

QF = uint32(80);
MEXstart = tic;
%% --------------------
% Run default DCTR extraction
% ---------------------
%fprintf('DCTR extraction - default settings');
F = DCTR({ImageSet}, QF);

%% --------------------
% Run custom PSRM extractions
% ---------------------
%{
config.T = uint32(3);    % number of histogram bins
config.q = single(1);    % bin width
fprintf('DCTR extraction - custom settings');
F = DCTR(ImageSet, QF, config);
%}

%% Results        
MEXend = toc(MEXstart);
%fprintf(' - DONE');
%fprintf('\n\nDCTR extracted %d images in %.2f seconds, in average %.2f seconds per image\n', numel(ImageSet), MEXend, MEXend / numel(ImageSet));