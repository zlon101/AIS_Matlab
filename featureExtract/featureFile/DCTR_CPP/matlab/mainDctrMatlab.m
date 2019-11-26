function F=mainDctrMatlab(imagePath)
%               EXAMPLE - USING "DCTR features"
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
%imagePath = fullfile('..', '..', 'image_dir', 'image.jpg');

% Specify JPEG image for feature extraction and its quality factor
quality_factor = 80;

%% ---------------------------
% DCTR extraction by MATLAB 
% ----------------------------
% fprintf('DCTR extraction');
I_STRUCT = jpeg_read(imagePath);

t_start = tic;
F = DCTR(I_STRUCT, quality_factor);
t_end = toc(t_start);

% fprintf(' - processed in %.2f seconds', t_end);

