function f = CFstar(IMAGE,QF)
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
% Contact: jan@kodovsky.com | fridrich@binghamton.edu | October 2011
%          http://dde.binghamton.edu/download/feature_extractors
% -------------------------------------------------------------------------
% Extracts CF* features [1] from the given IMAGE. Conveniently, all the
% features are extracted in the same order as they are listed in Table I in
% the publication [1]. The output is stored in a structured variable 'f'.
% The first half of the features are non-calibrated features as listed
% in Table I of [1]. The second half are the reference values that are
% extracted from the reference image (decompressed, cropped by 4x4,
% re-compressed with JPEG quality factor QF). In [1], we used QF=75 as both
% the original and the reference quality factor.  Total dimensionality of
% the features: 7850.
% -------------------------------------------------------------------------
% Input:  IMAGE .. path to the JPEG image
%         QF ..... JPEG quality factor of the reference image (we
%                   recommend to keep it the same as the original image QF)
% Output: f ...... extracted CF* features
% -------------------------------------------------------------------------
% IMPORTANT: For accessing DCT coefficients of JPEG images, we use Phil
% Sallee's Matlab JPEG toolbox (function jpeg_read) available at
% http://philsallee.com/jpegtbx/
% -------------------------------------------------------------------------
% Fintra = {Fh,Fd,Foh,Fx,Fod,Fm} - dim. 2375
% Finter = {Fih,Fid,Fis,Fix} - dim. 1550
% Fstar  = {Finter,Fintra} - dim. 3925
% CFstar = {Fstar} + reference version - dim. 7850
% -------------------------------------------------------------------------
% [1] Ensemble Classifiers for Steganalysis of Digital Media, J. Kodovsky,
% J. Fridrich and V. Holub, IEEE Transactions on Information Forensics and
% Security, 2011.
% -------------------------------------------------------------------------

if ~exist('QF','var'), QF = 75; end
X = DCTPlane(IMAGE);
T = 3;

% first half of the features (non-calibrated version)
f.Fh  = extract_submodels(X,1:20,T);
f.Fd  = extract_submodels(X,21:40,T);
f.Foh = extract_submodels(X,41:54,T);
f.Fx  = extract_submodels(X,55:63,T);
f.Fod = extract_submodels(X,114:130,T);
f.Fm  = extract_submodels(X,131:145,T);
f.Fih = extract_submodels(X,64:83,T);
f.Fid = extract_submodels(X,84:94,T);
f.Fis = extract_submodels(X,[95:99 101:104 106 107],T);
f.Fix = extract_submodels(X,146:165,T);

% reference version of the features (for calibration)
X = DCTPlane_reference(IMAGE,QF);
f.Fh_ref  = extract_submodels(X,1:20,T);
f.Fd_ref  = extract_submodels(X,21:40,T);
f.Foh_ref = extract_submodels(X,41:54,T);
f.Fx_ref  = extract_submodels(X,55:63,T);
f.Fod_ref = extract_submodels(X,114:130,T);
f.Fm_ref  = extract_submodels(X,131:145,T);
f.Fih_ref = extract_submodels(X,64:83,T);
f.Fid_ref = extract_submodels(X,84:94,T);
f.Fis_ref = extract_submodels(X,[95:99 101:104 106 107],T);
f.Fix_ref = extract_submodels(X,146:165,T);

function F = extract_submodels(X,IDs,T)

F = [];
for ID=IDs
    % interpret cooccurence ID as a list of DCT modes of interest
    target = interpret_ID(ID);
    % extract the cooccurence features
    columns = ExtractCoocColumns(X,target);
    f1 = extractCooccurencesFromColumns(columns,T);
    % extract 8x8 diagonally symmetric cooccurence features
    columns = ExtractCoocColumns(X,target(:,[2,1]));
    f2 = extractCooccurencesFromColumns(columns,T);
    % sign symmetrize
    f = sign_symmetrize_with_normalization(f1+f2);
    
    F = [F;f]; %#ok<AGROW>
end
function Plane=DCTPlane(path)
% loads DCT Plane of the given JPEG image + Quantization table
jobj=jpeg_read(path);
Plane=jobj.coef_arrays{1};
function Plane=DCTPlane_reference(path,QF)
I = imread(path);      % decompress into spatial domain
I = I(5:end-4,5:end-4); % crop by 4x4 pixels

TMP = ['img_' num2str(round(rand()*1e7)) num2str(round(rand()*1e7)) '.jpg'];
while exist(TMP,'file'), TMP = ['img_' num2str(round(rand()*1e7)) num2str(round(rand()*1e7)) '.jpg']; end
imwrite(I,TMP,'Quality',QF); % save as temporary jpeg image using imwrite
Plane = DCTPlane(TMP); % load DCT plane of the reference image
delete(TMP); % delete the temporary reference image
function target = interpret_ID(ID)
% takes cooccurence ID and returns list of DCT modes of interest

% ID 1-20, horizontally neighbouring pairs -
starting_points=[0 1 2 3 4 5;6 7 8 9 10 0;11 12 13 14 0 0;15 16 17 0 0 0;18 19 0 0 0 0;20 0 0 0 0 0];
[a,b] = ind2sub(size(starting_points),find(starting_points==ID));
if ~isempty(a), target = [a,b;a,b+1]; return; end
% ID 21-31, diagonally neighbouring pairs \
starting_points=[0 21 22 23 24 25;0 26 27 28 29 0;0 0 30 31 0 0];
[a,b] = ind2sub(size(starting_points),find(starting_points==(ID)));
if ~isempty(a), target = [a,b;a+1,b+1]; return; end
% ID 32-40, semidiagonally neighbouring pairs /
starting_points=[0 0 0 0 0;32 33 34 35 36;0 37 38 39 0;0 0 40 0 0];
[a,b] = ind2sub(size(starting_points),find(starting_points==(ID)));
if ~isempty(a), target = [a,b;a-1,b+1]; return; end
% ID 41-54, horizontally ob jedno
starting_points=[0 41 42 43 44;45 46 47 48 0;49 50 51 0 0;52 53 0 0 0;54 0 0 0 0];
[a,b] = ind2sub(size(starting_points),find(starting_points==(ID)));
if ~isempty(a), target = [a,b;a,b+2]; return; end
% ID 55-63 symmetric wrt diagonal
starting_points=[0 55 56 57 58 59;0 0 60 61 62 0;0 0 0 63 0 0];
[a,b] = ind2sub(size(starting_points),find(starting_points==(ID)));
if ~isempty(a), target = [a,b;b,a]; return; end
% ID 64-83 , horizontally neighbouring inter-block -
starting_points=[0 64 65 66 67 68;69 70 71 72 73 0;74 75 76 77 0 0;78 79 80 0 0 0;81 82 0 0 0 0;83 0 0 0 0 0];
[a,b] = ind2sub(size(starting_points),find(starting_points==ID));
if ~isempty(a), target = [a,b;a,b+8]; return; end
% ID 84-94 , inter-block diagonally
starting_points=[0 84 85 86 87 88;0 89 90 91 92 0;0 0 93 94 0 0];
[a,b] = ind2sub(size(starting_points),find(starting_points==ID));
if ~isempty(a), target = [a,b;a+8,b+8]; return; end
% ID 95-113,166 , inter-block semidiagonally
% duplicates: 100,166,105,108:113
% valid only: 95:99,101:104,106,107
starting_points=[0 95 96 97 98 99;100 101 102 103 104 0;166 105 106 107 0 0;108 109 110 0 0 0;111 112 0 0 0 0;113 0 0 0 0 0];
[a,b] = ind2sub(size(starting_points),find(starting_points==ID));
if ~isempty(a), target = [a+8,b;a,b+8]; return; end
% ID 114-121 , intra, diagonally ob jedno
starting_points=[0 114 115 116 117;0 118 119 120 0;0 0 121 0 0];
[a,b] = ind2sub(size(starting_points),find(starting_points==ID));
if ~isempty(a), target = [a,b;a+2,b+2]; return; end
% ID 122-130 , intra, semidiagonally ob jedno
starting_points=[0 0 0 0 0;0 0 0 0 0;122 123 124 125 126;0 127 128 129 0;0 0 130 0 0];
[a,b] = ind2sub(size(starting_points),find(starting_points==ID));
if ~isempty(a), target = [a,b;a-2,b+2]; return; end
% ID 131-145 , intra, horse-move
starting_points=[0 0 0 0 0;131 132 133 134 135;136 137 138 139 0;140 141 142 0 0;143 144 0 0 0;145 0 0 0 0];
[a,b] = ind2sub(size(starting_points),find(starting_points==ID));
if ~isempty(a), target = [a,b;a-1,b+2]; return; end
% ID 146-165 , inter, neighbouring diag. symmetric
starting_points=[0 146 147 148 149 150;151 152 153 154 155 0;156 157 158 159 0 0;160 161 162 0 0 0;163 164 0 0 0 0;165 0 0 0 0 0];
[a,b] = ind2sub(size(starting_points),find(starting_points==ID));
if ~isempty(a), target = [a,b;b,a+8]; return; end
function columns = ExtractCoocColumns(A,target)
% Take the target DCT modes and extracts their corresponding n-tuples from
% the DCT plane A. Store them as individual columns.
mask = getMask(target);
v = floor(size(A,1)/8)+1-(size(mask,1)/8); % number of vertical block shifts
h = floor(size(A,2)/8)+1-(size(mask,2)/8); % number of horizontal block shifts

for i=1:size(target,1)
    C = A(target(i,1)+(1:8:8*v)-1,target(i,2)+(1:8:8*h)-1);
    if ~exist('columns','var'),columns = zeros(numel(C),size(target,1)); end
    columns(:,i) = C(:);
end
function F = extractCooccurencesFromColumns(blocks,t)
% blocks = columns of values from which we want to extract the
% cooccurences. marginalize to [-t..t]. no normalization involved

order = size(blocks,2); % cooccurence order
blocks(blocks<-t) = -t; % marginalization
blocks(blocks>t) = t;   % marginalization

switch order
    case 1
        % 1st order cooccurence (histogram)
        F = zeros(2*t+1,1);
        % for loop is faster than hist function
        for i=-t:t
            F(i+t+1) = sum(blocks==i);
        end
    case 2
        % 2nd order cooccurence
        F = zeros(2*t+1,2*t+1);
        for i=-t:t
            fB = blocks(blocks(:,1)==i,2);
            if ~isempty(fB)
                for j=-t:t
                    F(i+t+1,j+t+1) = sum(fB==j);
                end
            end
        end
    case 3
        % 3rd order cooccurence
        F = zeros(2*t+1,2*t+1,2*t+1);
        for i=-t:t
            fB = blocks(blocks(:,1)==i,2:end);
            if ~isempty(fB)
                for j=-t:t
                    fB2 = fB(fB(:,1)==j,2:end);
                    if ~isempty(fB2)
                        for k=-t:t
                            F(i+t+1,j+t+1,k+t+1) = sum(fB2==k);
                        end
                    end
                end
            end
        end
    case 4
        % 4th order cooccurence
        F = zeros(2*t+1,2*t+1,2*t+1,2*t+1);
        for i=-t:t
            fB = blocks(blocks(:,1)==i,2:end);
            if ~isempty(fB)
                for j=-t:t
                    fB2 = fB(fB(:,1)==j,2:end);
                    if ~isempty(fB2)
                        for k=-t:t
                            fB3 = fB2(fB2(:,1)==k,2:end);
                            if ~isempty(fB3)
                                for l=-t:t
                                    F(i+t+1,j+t+1,k+t+1,l+t+1) = sum(fB3==l);
                                end
                            end
                        end
                    end
                end
            end
        end
end % switch
function mask = getMask(target)
% transform list of DCT modes of interest into a mask with all zeros and
% ones at the positions of those DCT modes of interest
x=8;y=8;
if sum(target(:,1)>8)>0 && sum(target(:,1)>16)==0, x=16; end
if sum(target(:,1)>16)>0 && sum(target(:,1)>24)==0, x=24; end
if sum(target(:,1)>24)>0 && sum(target(:,1)>32)==0, x=32; end
if sum(target(:,2)>8)>0 && sum(target(:,2)>16)==0, y=16; end
if sum(target(:,2)>16)>0 && sum(target(:,2)>24)==0, y=24; end
if sum(target(:,2)>24)>0 && sum(target(:,2)>32)==0, y=32; end

mask = zeros(x,y);
for i=1:size(target,1)
    mask(target(i,1),target(i,2)) = 1;
end
function f = sign_symmetrize_with_normalization(f)
% symmetrize cooccurence matrix by sign [x,y,z] = [-x,-y,-z] and output as
% a 1D vector (final form of the feature vector)
dim = numel(size(f));
sD = size(f,1);


% normalization at the beginning
f = f/sum(f(:));

switch dim
    case 1
        % 1d cooccurence (histogram) symmetrization
        center = (length(f)+1)/2;
        ID1 = 1:center-1;
        ID2 = length(f):-1:center+1;
    case 2
        % 2d cooccurence sign symmetrization
        REF = reshape(1:numel(f),sD,sD);
        DONE = zeros(size(REF));
        [ID1,ID2] = deal(zeros((numel(f)-1)/2,1));
        id = 0;
        for i=1:sD
            for j=1:sD
                if ~isequal([i,j],[sD-i+1,sD-j+1])&&~DONE(i,j)
                    id = id+1;
                    ID1(id) = REF(i,j);
                    ID2(id) = REF(sD-i+1,sD-j+1);
                    DONE(i,j)=1;
                    DONE(sD-i+1,sD-j+1)=1;
                end
            end
        end
    case 3
        % 3d cooccurence sign symmetrization
        REF = reshape(1:numel(f),sD,sD,sD);
        DONE = zeros(size(REF));
        [ID1,ID2] = deal(zeros((numel(f)-1)/2,1));
        id = 0;
        for i=1:sD
            for j=1:sD
                for k=1:sD
                    if ~isequal([i,j,k],[sD-i+1,sD-j+1,sD-k+1])&&~DONE(i,j,k)
                        id = id+1;
                        ID1(id) = REF(i,j,k);
                        ID2(id) = REF(sD-i+1,sD-j+1,sD-k+1);
                        DONE(i,j,k)=1;
                        DONE(sD-i+1,sD-j+1,sD-k+1)=1;
                    end
                end
            end
        end
    case 4
        % 4d cooccurence sign symmetrization
        REF = reshape(1:numel(f),sD,sD,sD,sD);
        DONE = zeros(size(REF));
        [ID1,ID2] = deal(zeros((numel(f)-1)/2,1));
        id = 0;
        for i=1:sD
            for j=1:sD
                for k=1:sD
                    for l=1:sD
                        if ~isequal([i,j,k,l],[sD-i+1,sD-j+1,sD-k+1,sD-l+1])&&~DONE(i,j,k,l)
                            id = id+1;
                            ID1(id) = REF(i,j,k,l);
                            ID2(id) = REF(sD-i+1,sD-j+1,sD-k+1,sD-l+1);
                            DONE(i,j,k,l)=1;
                            DONE(sD-i+1,sD-j+1,sD-k+1,sD-l+1)=1;
                        end
                    end
                end
            end
        end
end

% symmetrize
f = f(:);
f(ID1) = f(ID1)+f(ID2);
f(ID2) = [];