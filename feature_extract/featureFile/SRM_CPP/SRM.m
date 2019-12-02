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
% Extracts all 106 submodels presented in [1] as part of a rich model for
% steganalysis of digital images. All features are calculated in the
% spatial domain and are stored in a structured variable 'f'. For more
% deatils about the individual submodels, please see the publication [1].
% Total dimensionality of all 106 submodels is 34671.
% -------------------------------------------------------------------------
% F = SRM(imageSet) extracts the features from every image in 'imageSet'
% for each submodel using default settings
%
% F = SRM(imageSet, config)
% same as above but allows to configure the settings
%
%--------------------------------------------------------------------------
% Input format:
%   imageSet            - cell array of strings     - cell array of paths to the images
%
%   config.T			- int32		- default 2		- residual threshold
%   config.order		- int32		- default 4		- co-occurrence order
%   config.merge_spams	- logical	- default true	- if true then spam features are merged
%	config.symm_sign	- logical	- default true	- if true then spam symmetry is used
%	config.symm_reverse	- logical	- default true	- if true then reverse symmetry is used
%	config.symm_minmax	- logical	- default true	- if true then minmax symmetry is used
%	config.eraseLSB		- logical	- default false	- if true then all LSB are erased from the image
%	config.parity		- logical	- default false	- if true then parity residual is applied
%
% Output format:
%   F - Structure containing one fields for every submodel (the field is
%       identified by the submodel's name). Every field is a matrix of I*D
%       of type 'single' - I is number of images in imageSet, D is
%       feature dimensionality of the specific submodel (169 for SPAM
%       features, 325 for MinMax features default settings).
% -------------------------------------------------------------------------
% [1] Rich Models for Steganalysis of Digital Images, J. Fridrich and J.
% Kodovsky, IEEE Transactions on Information Forensics and Security, 2011.
% Under review. 
% http://dde.binghamton.edu/kodovsky/pdf/TIFS2012-SRM.pdf
% -------------------------------------------------------------------------
