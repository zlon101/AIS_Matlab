function f = PSRM3(IMAGE)
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
% Contact: vojtech_holub@yahoo.com | fridrich@binghamton.edu | January 2013
%          http://dde.binghamton.edu/download/feature_extractors
% -------------------------------------------------------------------------
% Extracts 39 submodels presented in [2] -- All features are
% calculated in the spatial domain and are stored in a structured variable
% 'f'. For more deatils about the individual submodels, please see the
% publication [1][2]. Total dimensionality of is 4290 x settings.projCount.
% -------------------------------------------------------------------------
% Input:  IMAGE ... path to the image or its pixel value matrix
% Output: f ...... extracted PSRM features in a structured format
% -------------------------------------------------------------------------
% [1] Random Projections of Residuals as an Alternative to Co-occurrences 
% in Steganalysis, V. Holub and J. Fridrich, Proc. SPIE, Electronic 
% Imaging, Media Watermarking, Security, and Forensics XV, 2013
%
% [2] Rich Models for Steganalysis of Digital Images, J. Fridrich and J.
% Kodovsky, IEEE Transactions on Information Forensics and Security, 2011.
% Under review.
% -------------------------------------------------------------------------

% Setting for PSRM3
settings.projCount = 3; 

% Setting for PSRM8
%settings.projCount = 8; 

settings.B = 2;
settings.qBins = 5;
settings.seedIndex = 1;

settings.neighborhoods{1} = [1, 1, 1, 1];
settings.neighborhoods{2} = [1, 1, 1, 1, 1, 1, 1, 1];
settings.neighborhoods{3} = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
settings.neighborhoods{4} = [1, 1, 0, 0; 1, 1, 1, 0; 0, 1, 1, 1; 0, 0, 1, 1];
settings.neighborhoods{5} = [1, 1; 1, 1];
settings.neighborhoods{6} = [1, 1, 1; 1, 1, 1; 1, 1, 1];
settings.neighborhoods{7} = [1, 1, 1, 1; 1, 1, 1, 1; 1, 1, 1, 1; 1, 1, 1, 1];
settings.neighborhoods{8} = [1, 1, 1, 1, 1; 1, 1, 1, 1, 1; 1, 1, 1, 1, 1; 1, 1, 1, 1, 1; 1, 1, 1, 1, 1];
settings.neighborhoods{9} = [1, 1, 1, 1; 1, 1, 1, 1];
settings.neighborhoods{10} = [1, 0, 0, 0; 1, 1, 0, 0; 1, 1, 1, 0; 1, 1, 1, 1];
settings.neighborhoods{11} = [0, 0, 1, 0, 0; 0, 0, 1, 0, 0; 1, 1, 1, 1, 1; 0, 0, 1, 0, 0; 0, 0, 1, 0, 0];

if ischar(IMAGE)  
    X = double(imread(IMAGE));
else
    X = double(IMAGE);
end    

f = post_processing(all1st(X,1),'f1');      % 1st order
f = post_processing(all2nd(X,2),'f2',f);    % 2nd order
f = post_processing(all3rd(X,3),'f3',f);    % 3rd order
f = post_processing(all3x3(X,4),'f3x3', f);  % 3x3
f = post_processing(all5x5(X,12),'f5x5', f); % 5x5

function RESULT = post_processing(DATA,f,RESULT)

Ss = fieldnames(DATA);
for Sid = 1:length(Ss)
    VARNAME = [f '_' Ss{Sid}];
    eval(['RESULT.' VARNAME ' = reshape(single(DATA.' Ss{Sid} '),1,[]);' ])
end

% symmetrize
L = fieldnames(RESULT);
for i=1:length(L)
    name = L{i}; % feature name
    if name(1)=='s', continue; end
    [T,N] = parse_feaname(name);
    if strcmp(T,''), continue; end
    % symmetrization
    if strcmp(N(1:3),'min') || strcmp(N(1:3),'max')
        % minmax symmetrization
        OUT = ['s' T(2:end) '_minmax' N(4:end)];
        if isfield(RESULT,OUT), continue; end
        Fmin = []; Fmax = [];
        eval(['Fmin = RESULT.' strrep(name,'max','min') ';']);
        eval(['Fmax = RESULT.' strrep(name,'min','max') ';']);
        F = mergeMinMax(Fmin, Fmax);
        eval(['RESULT.' OUT ' = single(F);' ]);
    elseif strcmp(N(1:4),'spam')
        % spam symmetrization
        OUT = ['s' T(2:end) '_' N];
        if isfield(RESULT,OUT), continue; end
        eval(['RESULT.' OUT ' = single(RESULT.' name ');' ]);
    end
end
% delete RESULT.f*
L = fieldnames(RESULT);
for i=1:length(L)
    name = L{i}; % feature name
    if name(1)=='f'
        RESULT = rmfield(RESULT,name);
    end
end
% merge spam features
L = fieldnames(RESULT);
for i=1:length(L)
    name = L{i}; % feature name
    [T,N] = parse_feaname(name);
    if ~strcmp(N(1:4),'spam'), continue; end
    if strcmp(T,''), continue; end
    if strcmp(N(end),'v')||(strcmp(N,'spam11')&&strcmp(T,'s5x5'))
    elseif strcmp(N(end),'h')
        % h+v union
        OUT = [T '_' N 'v'];
        if isfield(RESULT,OUT), continue; end
        name2 = strrep(name,'h','v');
        Fh = []; Fv = [];
        eval(['Fh = RESULT.' name ';']);
        eval(['Fv = RESULT.' name2 ';']);
        eval(['RESULT.' OUT ' = [Fh Fv];']);
        RESULT = rmfield(RESULT,name);
        RESULT = rmfield(RESULT,name2);
    elseif strcmp(N,'spam11')
        % KBKV creation
        OUT = ['s35_' N];
        if isfield(RESULT,OUT), continue; end
        name1 = strrep(name,'5x5','3x3');
        name2 = strrep(name,'3x3','5x5');
        if ~isfield(RESULT,name1), continue; end
        if ~isfield(RESULT,name2), continue; end
        F_KB = []; F_KV = [];
        eval(['F_KB = RESULT.' name1 ';']);
        eval(['F_KV = RESULT.' name2 ';']);
        eval(['RESULT.' OUT ' = [F_KB F_KV];']);
        RESULT = rmfield(RESULT,name1);
        RESULT = rmfield(RESULT,name2);
    end
end

end

function [T,N] = parse_feaname(name)
[T,N] = deal('');
S = strfind(name,'_'); if length(S)~=1, return; end
T = name(1:S-1);
N = name(S+1:end);

end

function g = all1st(X,q)
%
% X must be a matrix of doubles or singles (the image) and q is the 
% quantization step (any positive number).
%
% Recommended values of q are c, 1.5c, 2c, where c is the central
% coefficient in the differential (at X(I,J)).
%
% This function outputs co-occurrences of ALL 1st-order residuals
% listed in Figure 1 in our journal HUGO paper (version from June 14), 
% including the naming convention.
%
% List of outputted features:
%
% 1a) spam14h
% 1b) spam14v (orthogonal-spam)
% 1c) minmax22v
% 1d) minmax24
% 1e) minmax34v
% 1f) minmax41
% 1g) minmax34
% 1h) minmax48h
% 1i) minmax54
%
% Naming convention:
%
% name = {type}{f}{sigma}{scan}
% type \in {spam, minmax}
% f \in {1,2,3,4,5} number of filters that are "minmaxed"
% sigma \in {1,2,3,4,8} symmetry index
% scan \in {h,v,\emptyset} scan of the cooc matrix (empty = sum of both 
% h and v scans).
%
% All odd residuals are implemented the same way simply by
% narrowing the range for I and J and replacing the residuals --
% -- they should "stick out" (trcet) in the same direction as 
% the 1st order ones. For example, for the 3rd order:
%
% RU = -X(I-2,J+2)+3*X(I-1,J+1)-3*X(I,J)+X(I+1,J-1); ... etc.
%
% Note1: The term X(I,J) should always have the "-" sign.
% Note2: This script does not include s, so, cout, cin versions (weak).
% This function calls Cooc.m and Quant.m
[M N] = size(X); [I,J] = deal(2:M-1,2:N-1);
% Variable names are self-explanatory (R = right, U = up, L = left, D = down)
[R,L,U,D]  = deal(X(I,J+1)-X(I,J),X(I,J-1)-X(I,J),X(I-1,J)-X(I,J),X(I+1,J)-X(I,J));
[RU,LU,RD,LD] = deal(X(I-1,J+1)-X(I,J),X(I-1,J-1)-X(I,J),X(I+1,J+1)-X(I,J),X(I+1,J-1)-X(I,J));
% minmax22h -- to be symmetrized as mnmx, directional, hv-nonsymmetrical.
[RLq_min,UDq_min,RLq_max,UDq_max] = deal(min(R,L),min(U,D),max(R,L),max(U,D));
g.min22h = reshape(ProjHistMinMax(RLq_min,'hor',q) + ProjHistMinMax(UDq_min,'ver',q),[],1);
g.max22h = reshape(ProjHistMinMax(RLq_max,'hor',q) + ProjHistMinMax(UDq_max,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax34h -- to be symmetrized as mnmx, directional, hv-nonsymmetrical
[Uq_min,Rq_min,Dq_min,Lq_min] = deal(min(min(L,U),R),min(min(U,R),D),min(min(R,D),L),min(min(D,L),U));
[Uq_max,Rq_max,Dq_max,Lq_max] = deal(max(max(L,U),R),max(max(U,R),D),max(max(R,D),L),max(max(D,L),U));
g.min34h = reshape(ProjHistMinMax(Uq_min,'hor',q) + ProjHistMinMax(Dq_min,'hor',q) + ProjHistMinMax(Lq_min,'ver',q) + ProjHistMinMax(Rq_min,'ver',q),[],1);
g.max34h = reshape(ProjHistMinMax(Uq_max,'hor',q) + ProjHistMinMax(Dq_max,'hor',q) + ProjHistMinMax(Lq_max,'ver',q) + ProjHistMinMax(Rq_max,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
% spam14h/v -- to be symmetrized as spam, directional, hv-nonsymmetrical
g.spam14h = reshape(ProjHistSpam(R,'hor',q) + ProjHistSpam(U,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
g.spam14v = reshape(ProjHistSpam(R,'ver',q) + ProjHistSpam(U,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax22v -- to be symmetrized as mnmx, directional, hv-nonsymmetrical. Good with higher-order residuals! Note: 22h is bad (too much neighborhood overlap).
g.min22v = reshape(ProjHistMinMax(RLq_min,'ver',q) + ProjHistMinMax(UDq_min,'hor',q),[],1);
g.max22v = reshape(ProjHistMinMax(RLq_max,'ver',q) + ProjHistMinMax(UDq_max,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax24 -- to be symmetrized as mnmx, directional, hv-symmetrical. Darn good, too.
[RUq_min,RDq_min,LUq_min,LDq_min] = deal(min(R,U),min(R,D),min(L,U),min(L,D));
[RUq_max,RDq_max,LUq_max,LDq_max] = deal(max(R,U),max(R,D),max(L,U),max(L,D));
g.min24 = reshape(ProjHistMinMax(RUq_min,'hor',q) + ProjHistMinMax(RDq_min,'hor',q) + ProjHistMinMax(LUq_min,'hor',q) + ProjHistMinMax(LDq_min,'hor',q) +...
                  ProjHistMinMax(RUq_min,'ver',q) + ProjHistMinMax(RDq_min,'ver',q) + ProjHistMinMax(LUq_min,'ver',q) + ProjHistMinMax(LDq_min,'ver',q),[],1);
g.max24 = reshape(ProjHistMinMax(RUq_max,'hor',q) + ProjHistMinMax(RDq_max,'hor',q) + ProjHistMinMax(LUq_max,'hor',q) + ProjHistMinMax(LDq_max,'hor',q) +...
                  ProjHistMinMax(RUq_max,'ver',q) + ProjHistMinMax(RDq_max,'ver',q) + ProjHistMinMax(LUq_max,'ver',q) + ProjHistMinMax(LDq_max,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax34v -- v works well, h does not, to be symmetrized as mnmx, directional, hv-nonsymmetrical
g.min34v = reshape(ProjHistMinMax(Uq_min,'ver',q) + ProjHistMinMax(Dq_min,'ver',q) + ProjHistMinMax(Rq_min,'hor',q) + ProjHistMinMax(Lq_min,'hor',q),[],1);
g.max34v = reshape(ProjHistMinMax(Uq_max,'ver',q) + ProjHistMinMax(Dq_max,'ver',q) + ProjHistMinMax(Rq_max,'hor',q) + ProjHistMinMax(Lq_max,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax41 -- to be symmetrized as mnmx, non-directional, hv-symmetrical
[R_min,R_max] = deal(min(RLq_min,UDq_min),max(RLq_max,UDq_max));
g.min41 = reshape(ProjHistMinMax(R_min,'hor',q) + ProjHistMinMax(R_min,'ver',q),[],1);
g.max41 = reshape(ProjHistMinMax(R_max,'hor',q) + ProjHistMinMax(R_max,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax34 -- good, to be symmetrized as mnmx, directional, hv-symmetrical
[RUq_min,RDq_min,LUq_min,LDq_min] = deal(min(RUq_min,RU),min(RDq_min,RD),min(LUq_min,LU),min(LDq_min,LD));
[RUq_max,RDq_max,LUq_max,LDq_max] = deal(max(RUq_max,RU),max(RDq_max,RD),max(LUq_max,LU),max(LDq_max,LD));
g.min34 = reshape(ProjHistMinMax(RUq_min,'hor',q) + ProjHistMinMax(RDq_min,'hor',q) + ProjHistMinMax(LUq_min,'hor',q) + ProjHistMinMax(LDq_min,'hor',q) + ...
                  ProjHistMinMax(RUq_min,'ver',q) + ProjHistMinMax(RDq_min,'ver',q) + ProjHistMinMax(LUq_min,'ver',q) + ProjHistMinMax(LDq_min,'ver',q),[],1);
g.max34 = reshape(ProjHistMinMax(RUq_max,'hor',q) + ProjHistMinMax(RDq_max,'hor',q) + ProjHistMinMax(LUq_max,'hor',q) + ProjHistMinMax(LDq_max,'hor',q) + ...
                  ProjHistMinMax(RUq_max,'ver',q) + ProjHistMinMax(RDq_max,'ver',q) + ProjHistMinMax(LUq_max,'ver',q) + ProjHistMinMax(LDq_max,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax48h -- h better than v, to be symmetrized as mnmx, directional, hv-nonsymmetrical. 48v is almost as good as 48h; for 3rd-order but weaker for 1st-order. Here, I am outputting both but Figure 1 in our paper lists only 48h.
[RUq_min2,RDq_min2,LDq_min2,LUq_min2] = deal(min(RUq_min,LU),min(RDq_min,RU),min(LDq_min,RD),min(LUq_min,LD));
[RUq_min3,RDq_min3,LDq_min3,LUq_min3] = deal(min(RUq_min,RD),min(RDq_min,LD),min(LDq_min,LU),min(LUq_min,RU));
g.min48h = reshape(ProjHistMinMax(RUq_min2,'hor',q) + ProjHistMinMax(LDq_min2,'hor',q) + ProjHistMinMax(RDq_min3,'hor',q) + ProjHistMinMax(LUq_min3,'hor',q) + ...
                   ProjHistMinMax(RDq_min2,'ver',q) + ProjHistMinMax(LUq_min2,'ver',q) + ProjHistMinMax(RUq_min3,'ver',q) + ProjHistMinMax(LDq_min3,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
g.min48v = reshape(ProjHistMinMax(RDq_min2,'hor',q) + ProjHistMinMax(LUq_min2,'hor',q) + ProjHistMinMax(RUq_min3,'hor',q) + ProjHistMinMax(LDq_min3,'hor',q) + ...
                   ProjHistMinMax(RUq_min2,'ver',q) + ProjHistMinMax(LDq_min2,'ver',q) + ProjHistMinMax(RDq_min3,'ver',q) + ProjHistMinMax(LUq_min3,'ver',q),[],1);settings.seedIndex=settings.seedIndex-1;
[RUq_max2,RDq_max2,LDq_max2,LUq_max2] = deal(max(RUq_max,LU),max(RDq_max,RU),max(LDq_max,RD),max(LUq_max,LD));
[RUq_max3,RDq_max3,LDq_max3,LUq_max3] = deal(max(RUq_max,RD),max(RDq_max,LD),max(LDq_max,LU),max(LUq_max,RU));
g.max48h = reshape(ProjHistMinMax(RUq_max2,'hor',q) + ProjHistMinMax(LDq_max2,'hor',q) + ProjHistMinMax(RDq_max3,'hor',q) + ProjHistMinMax(LUq_max3,'hor',q) + ...
                   ProjHistMinMax(RDq_max2,'ver',q) + ProjHistMinMax(LUq_max2,'ver',q) + ProjHistMinMax(RUq_max3,'ver',q) + ProjHistMinMax(LDq_max3,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
g.max48v = reshape(ProjHistMinMax(RDq_max2,'hor',q) + ProjHistMinMax(LUq_max2,'hor',q) + ProjHistMinMax(RUq_max3,'hor',q) + ProjHistMinMax(LDq_max3,'hor',q) + ...
                   ProjHistMinMax(RUq_max2,'ver',q) + ProjHistMinMax(LDq_max2,'ver',q) + ProjHistMinMax(RDq_max3,'ver',q) + ProjHistMinMax(LUq_max3,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax54 -- to be symmetrized as mnmx, directional, hv-symmetrical
[RUq_min4,RDq_min4,LDq_min4,LUq_min4] = deal(min(RUq_min2,RD),min(RDq_min2,LD),min(LDq_min2,LU),min(LUq_min2,RU));
[RUq_min5,RDq_min5,LDq_min5,LUq_min5] = deal(min(RUq_min3,LU),min(RDq_min3,RU),min(LDq_min3,RD),min(LUq_min3,LD));
g.min54 = reshape(ProjHistMinMax(RUq_min4,'hor',q) + ProjHistMinMax(LDq_min4,'hor',q) + ProjHistMinMax(RDq_min5,'hor',q) + ProjHistMinMax(LUq_min5,'hor',q) +  ...
                  ProjHistMinMax(RDq_min4,'ver',q) + ProjHistMinMax(LUq_min4,'ver',q) + ProjHistMinMax(RUq_min5,'ver',q) + ProjHistMinMax(LDq_min5,'ver',q),[],1);
[RUq_max4,RDq_max4,LDq_max4,LUq_max4] = deal(max(RUq_max2,RD),max(RDq_max2,LD),max(LDq_max2,LU),max(LUq_max2,RU));
[RUq_max5,RDq_max5,LDq_max5,LUq_max5] = deal(max(RUq_max3,LU),max(RDq_max3,RU),max(LDq_max3,RD),max(LUq_max3,LD));
g.max54 = reshape(ProjHistMinMax(RUq_max4,'hor',q) + ProjHistMinMax(LDq_max4,'hor',q) + ProjHistMinMax(RDq_max5,'hor',q) + ProjHistMinMax(LUq_max5,'hor',q) +  ...
                  ProjHistMinMax(RDq_max4,'ver',q) + ProjHistMinMax(LUq_max4,'ver',q) + ProjHistMinMax(RUq_max5,'ver',q) + ProjHistMinMax(LDq_max5,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;

end

function g = all2nd(X,q)
%
% X must be a matrix of doubles or singles (the image) and q is the 
% quantization step (any positive number).
%
% Recommended values of q are c, 1.5c, 2c, where c is the central
% coefficient in the differential (at X(I,J)).
%
% This function outputs co-occurrences of ALL 2nd-order residuals
% listed in Figure 1 in our journal HUGO paper (version from June 14), 
% including the naming convention.
%
% List of outputted features:
%
% 1a) spam12h
% 1b) spam12v (orthogonal-spam)
% 1c) minmax21
% 1d) minmax41
% 1e) minmax24h (24v is also outputted but not listed in Figure 1)
% 1f) minmax32
%
% Naming convention:
%
% name = {type}{f}{sigma}{scan}
% type \in {spam, minmax}
% f \in {1,2,3,4,5} number of filters that are "minmaxed"
% sigma \in {1,2,3,4,8} symmetry index
% scan \in {h,v,\emptyset} scan of the cooc matrix (empty = sum of both 
% h and v scans).
%
% All even residuals are implemented the same way simply by
% narrowing the range for I and J and replacing the residuals.
%
% Note1: The term X(I,J) should always have the "-" sign.
% Note2: This script does not include s, so, cout, cin versions (weak).
%
% This function calls Residual.m, Cooc.m, and Quant.m

% 2nd-order residuals are implemented using Residual.m
[Dh,Dv,Dd,Dm] = deal(Residual(X,2,'hor'),Residual(X,2,'ver'),Residual(X,2,'diag'),Residual(X,2,'mdiag'));
% spam12h/v
g.spam12h = reshape(ProjHistSpam(Dh,'hor',q) + ProjHistSpam(Dv,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
g.spam12v = reshape(ProjHistSpam(Dh,'ver',q) + ProjHistSpam(Dv,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax21
[Dmin,Dmax] = deal(min(Dh,Dv),max(Dh,Dv));
g.min21 = reshape(ProjHistMinMax(Dmin,'hor',q) + ProjHistMinMax(Dmin,'ver',q),[],1);
g.max21 = reshape(ProjHistMinMax(Dmax,'hor',q) + ProjHistMinMax(Dmax,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax41   
[Dmin2,Dmax2] = deal(min(Dmin,min(Dd,Dm)),max(Dmax,max(Dd,Dm)));
g.min41 = reshape(ProjHistMinMax(Dmin2,'hor',q) + ProjHistMinMax(Dmin2,'ver',q),[],1);
g.max41 = reshape(ProjHistMinMax(Dmax2,'hor',q) + ProjHistMinMax(Dmax2,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax32 -- good, directional, hv-symmetrical, to be symmetrized as mnmx
[RUq_min,RDq_min] = deal(min(Dmin,Dm),min(Dmin,Dd));
[RUq_max,RDq_max] = deal(max(Dmax,Dm),max(Dmax,Dd));
g.min32 = reshape(ProjHistMinMax(RUq_min,'hor',q) + ProjHistMinMax(RDq_min,'hor',q) + ProjHistMinMax(RUq_min,'ver',q) + ProjHistMinMax(RDq_min,'ver',q),[],1);
g.max32 = reshape(ProjHistMinMax(RUq_max,'hor',q) + ProjHistMinMax(RDq_max,'hor',q) + ProjHistMinMax(RUq_max,'ver',q) + ProjHistMinMax(RDq_max,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax24h,v -- both "not bad," h slightly better, directional, hv-nonsymmetrical, to be symmetrized as mnmx
[RUq_min2,RDq_min2,RUq_min3,LUq_min3] = deal(min(Dm,Dh),min(Dd,Dh),min(Dm,Dv),min(Dd,Dv));
g.min24h = reshape(ProjHistMinMax(RUq_min2,'hor',q) + ProjHistMinMax(RDq_min2,'hor',q) + ProjHistMinMax(RUq_min3,'ver',q) + ProjHistMinMax(LUq_min3,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
g.min24v = reshape(ProjHistMinMax(RUq_min2,'ver',q) + ProjHistMinMax(RDq_min2,'ver',q) + ProjHistMinMax(RUq_min3,'hor',q) + ProjHistMinMax(LUq_min3,'hor',q),[],1);settings.seedIndex=settings.seedIndex-1;
[RUq_max2,RDq_max2,RUq_max3,LUq_max3] = deal(max(Dm,Dh),max(Dd,Dh),max(Dm,Dv),max(Dd,Dv));
g.max24h = reshape(ProjHistMinMax(RUq_max2,'hor',q) + ProjHistMinMax(RDq_max2,'hor',q) + ProjHistMinMax(RUq_max3,'ver',q) + ProjHistMinMax(LUq_max3,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
g.max24v = reshape(ProjHistMinMax(RUq_max2,'ver',q) + ProjHistMinMax(RDq_max2,'ver',q) + ProjHistMinMax(RUq_max3,'hor',q) + ProjHistMinMax(LUq_max3,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;

end

function g = all3rd(X,q)
%
% X must be a matrix of doubles or singles (the image) and q is the 
% quantization step (any positive number).
%
% Recommended values of q are c, 1.5c, 2c, where c is the central
% coefficient in the differential (at X(I,J)).
%
% This function outputs co-occurrences of ALL 3rd-order residuals
% listed in Figure 1 in our journal HUGO paper (version from June 14), 
% including the naming convention.
%
% List of outputted features:
%
% 1a) spam14h
% 1b) spam14v (orthogonal-spam)
% 1c) minmax22v
% 1d) minmax24
% 1e) minmax34v
% 1f) minmax41
% 1g) minmax34
% 1h) minmax48h
% 1i) minmax54
%
% Naming convention:
%
% name = {type}{f}{sigma}{scan}
% type \in {spam, minmax}
% f \in {1,2,3,4,5} number of filters that are "minmaxed"
% sigma \in {1,2,3,4,8} symmetry index
% scan \in {h,v,\emptyset} scan of the cooc matrix (empty = sum of both 
% h and v scans).
%
% All odd residuals are implemented the same way simply by
% narrowing the range for I and J and replacing the residuals --
% -- they should "stick out" (trcet) in the same direction as 
% the 1st order ones. For example, for the 3rd order:
%
% RU = -X(I-2,J+2)+3*X(I-1,J+1)-3*X(I,J)+X(I+1,J-1); ... etc.
%
% Note1: The term X(I,J) should always have the "-" sign.
% Note2: This script does not include s, so, cout, cin versions (weak).

[M N] = size(X); [I,J,] = deal(3:M-2,3:N-2);
[R,L,U,D] = deal(-X(I,J+2)+3*X(I,J+1)-3*X(I,J)+X(I,J-1),-X(I,J-2)+3*X(I,J-1)-3*X(I,J)+X(I,J+1),-X(I-2,J)+3*X(I-1,J)-3*X(I,J)+X(I+1,J),-X(I+2,J)+3*X(I+1,J)-3*X(I,J)+X(I-1,J));
[RU,LU,RD,LD] = deal(-X(I-2,J+2)+3*X(I-1,J+1)-3*X(I,J)+X(I+1,J-1),-X(I-2,J-2)+3*X(I-1,J-1)-3*X(I,J)+X(I+1,J+1),-X(I+2,J+2)+3*X(I+1,J+1)-3*X(I,J)+X(I-1,J-1),-X(I+2,J-2)+3*X(I+1,J-1)-3*X(I,J)+X(I-1,J+1));

% minmax22h -- to be symmetrized as mnmx, directional, hv-nonsymmetrical
[RLq_min,UDq_min,RLq_max,UDq_max] = deal(min(R,L),min(U,D),max(R,L),max(U,D));
g.min22h = reshape(ProjHistMinMax(RLq_min,'hor',q) + ProjHistMinMax(UDq_min,'ver',q),[],1);
g.max22h = reshape(ProjHistMinMax(RLq_max,'hor',q) + ProjHistMinMax(UDq_max,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax34h -- to be symmetrized as mnmx, directional, hv-nonsymmetrical
[Uq_min,Rq_min,Dq_min,Lq_min] = deal(min(RLq_min,U),min(UDq_min,R),min(RLq_min,D),min(UDq_min,L));
[Uq_max,Rq_max,Dq_max,Lq_max] = deal(max(RLq_max,U),max(UDq_max,R),max(RLq_max,D),max(UDq_max,L));
g.min34h = reshape(ProjHistMinMax(Uq_min,'hor',q) + ProjHistMinMax(Dq_min,'hor',q) + ProjHistMinMax(Lq_min,'ver',q) + ProjHistMinMax(Rq_min,'ver',q),[],1);
g.max34h = reshape(ProjHistMinMax(Uq_max,'hor',q) + ProjHistMinMax(Dq_max,'hor',q) + ProjHistMinMax(Lq_max,'ver',q) + ProjHistMinMax(Rq_max,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
% spam14h,v -- to be symmetrized as spam, directional, hv-nonsymmetrical
g.spam14h = reshape(ProjHistSpam(R,'hor',q) + ProjHistSpam(U,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
g.spam14v = reshape(ProjHistSpam(R,'ver',q) + ProjHistSpam(U,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax22v -- to be symmetrized as mnmx, directional, hv-nonsymmetrical. Good with higher-order residuals! Note: 22h is bad (too much neighborhood overlap).
g.min22v = reshape(ProjHistMinMax(RLq_min,'ver',q) + ProjHistMinMax(UDq_min,'hor',q),[],1);
g.max22v = reshape(ProjHistMinMax(RLq_max,'ver',q) + ProjHistMinMax(UDq_max,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax24 -- to be symmetrized as mnmx, directional, hv-symmetrical Note: Darn good, too.
[RUq_min,RDq_min,LUq_min,LDq_min] = deal(min(R,U),min(R,D),min(L,U),min(L,D));
[RUq_max,RDq_max,LUq_max,LDq_max] = deal(max(R,U),max(R,D),max(L,U),max(L,D));
g.min24 = reshape(ProjHistMinMax(RUq_min,'hor',q) + ProjHistMinMax(RDq_min,'hor',q) + ProjHistMinMax(LUq_min,'hor',q) + ProjHistMinMax(LDq_min,'hor',q) +...
                  ProjHistMinMax(RUq_min,'ver',q) + ProjHistMinMax(RDq_min,'ver',q) + ProjHistMinMax(LUq_min,'ver',q) + ProjHistMinMax(LDq_min,'ver',q),[],1);
g.max24 = reshape(ProjHistMinMax(RUq_max,'hor',q) + ProjHistMinMax(RDq_max,'hor',q) + ProjHistMinMax(LUq_max,'hor',q) + ProjHistMinMax(LDq_max,'hor',q) +...
                  ProjHistMinMax(RUq_max,'ver',q) + ProjHistMinMax(RDq_max,'ver',q) + ProjHistMinMax(LUq_max,'ver',q) + ProjHistMinMax(LDq_max,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax34v -- v works well, h does not, to be symmetrized as mnmx, directional, hv-nonsymmetrical
g.min34v = reshape(ProjHistMinMax(Uq_min,'ver',q) + ProjHistMinMax(Dq_min,'ver',q) + ProjHistMinMax(Rq_min,'hor',q) + ProjHistMinMax(Lq_min,'hor',q),[],1);
g.max34v = reshape(ProjHistMinMax(Uq_max,'ver',q) + ProjHistMinMax(Dq_max,'ver',q) + ProjHistMinMax(Rq_max,'hor',q) + ProjHistMinMax(Lq_max,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax41 -- unknown performance as of 6/14/11, to be symmetrized as mnmx, non-directional, hv-symmetrical
[R_min,R_max] = deal(min(RUq_min,LDq_min),max(RUq_max,LDq_max));
g.min41 = reshape(ProjHistMinMax(R_min,'hor',q) + ProjHistMinMax(R_min,'ver',q),[],1);
g.max41 = reshape(ProjHistMinMax(R_max,'hor',q) + ProjHistMinMax(R_max,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax34 -- good, to be symmetrized as mnmx, directional, hv-symmetrical
[RUq_min2,RDq_min2,LUq_min2,LDq_min2] = deal(min(RUq_min,RU),min(RDq_min,RD),min(LUq_min,LU),min(LDq_min,LD));
[RUq_max2,RDq_max2,LUq_max2,LDq_max2] = deal(max(RUq_max,RU),max(RDq_max,RD),max(LUq_max,LU),max(LDq_max,LD));
g.min34 = reshape(ProjHistMinMax(RUq_min2,'hor',q) + ProjHistMinMax(RDq_min2,'hor',q) + ProjHistMinMax(LUq_min2,'hor',q) + ProjHistMinMax(LDq_min2,'hor',q) + ...
                  ProjHistMinMax(RUq_min2,'ver',q) + ProjHistMinMax(RDq_min2,'ver',q) + ProjHistMinMax(LUq_min2,'ver',q) + ProjHistMinMax(LDq_min2,'ver',q),[],1);
g.max34 = reshape(ProjHistMinMax(RUq_max2,'hor',q) + ProjHistMinMax(RDq_max2,'hor',q) + ProjHistMinMax(LUq_max2,'hor',q) + ProjHistMinMax(LDq_max2,'hor',q) + ...
                  ProjHistMinMax(RUq_max2,'ver',q) + ProjHistMinMax(RDq_max2,'ver',q) + ProjHistMinMax(LUq_max2,'ver',q) + ProjHistMinMax(LDq_max2,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax48h -- h better than v, to be symmetrized as mnmx, directional, hv-nonsymmetrical. 48v is almost as good as 48h for 3rd-order but weaker for 1st-order. Here, I am outputting both but Figure 1 in our paper lists only 48h.
[RUq_min3,RDq_min3,LDq_min3,LUq_min3] = deal(min(RUq_min2,LU),min(RDq_min2,RU),min(LDq_min2,RD),min(LUq_min2,LD));
[RUq_min4,RDq_min4,LDq_min4,LUq_min4] = deal(min(RUq_min2,RD),min(RDq_min2,LD),min(LDq_min2,LU),min(LUq_min2,RU));
g.min48h = reshape(ProjHistMinMax(RUq_min3,'hor',q) + ProjHistMinMax(LDq_min3,'hor',q) + ProjHistMinMax(RDq_min4,'hor',q) + ProjHistMinMax(LUq_min4,'hor',q) + ...
                   ProjHistMinMax(RDq_min3,'ver',q) + ProjHistMinMax(LUq_min3,'ver',q) + ProjHistMinMax(RUq_min4,'ver',q) + ProjHistMinMax(LDq_min4,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
g.min48v = reshape(ProjHistMinMax(RUq_min3,'ver',q) + ProjHistMinMax(LDq_min3,'ver',q) + ProjHistMinMax(RDq_min4,'ver',q) + ProjHistMinMax(LUq_min4,'ver',q) + ...
                   ProjHistMinMax(RDq_min3,'hor',q) + ProjHistMinMax(LUq_min3,'hor',q) + ProjHistMinMax(RUq_min4,'hor',q) + ProjHistMinMax(LDq_min4,'hor',q),[],1);settings.seedIndex=settings.seedIndex-1;
[RUq_max3,RDq_max3,LDq_max3,LUq_max3] = deal(max(RUq_max2,LU),max(RDq_max2,RU),max(LDq_max2,RD),max(LUq_max2,LD));
[RUq_max4,RDq_max4,LDq_max4,LUq_max4] = deal(max(RUq_max2,RD),max(RDq_max2,LD),max(LDq_max2,LU),max(LUq_max2,RU));
g.max48h = reshape(ProjHistMinMax(RUq_max3,'hor',q) + ProjHistMinMax(LDq_max3,'hor',q) + ProjHistMinMax(RDq_max4,'hor',q) + ProjHistMinMax(LUq_max4,'hor',q) + ...
                   ProjHistMinMax(RDq_max3,'ver',q) + ProjHistMinMax(LUq_max3,'ver',q) + ProjHistMinMax(RUq_max4,'ver',q) + ProjHistMinMax(LDq_max4,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
g.max48v = reshape(ProjHistMinMax(RUq_max3,'ver',q) + ProjHistMinMax(LDq_max3,'ver',q) + ProjHistMinMax(RDq_max4,'ver',q) + ProjHistMinMax(LUq_max4,'ver',q) + ...
                   ProjHistMinMax(RDq_max3,'hor',q) + ProjHistMinMax(LUq_max3,'hor',q) + ProjHistMinMax(RUq_max4,'hor',q) + ProjHistMinMax(LDq_max4,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax54 -- to be symmetrized as mnmx, directional, hv-symmetrical
[RUq_min5,RDq_min5,LDq_min5,LUq_min5] = deal(min(RUq_min3,RD),min(RDq_min3,LD),min(LDq_min3,LU),min(LUq_min3,RU));
[RUq_max5,RDq_max5,LDq_max5,LUq_max5] = deal(max(RUq_max3,RD),max(RDq_max3,LD),max(LDq_max3,LU),max(LUq_max3,RU));
g.min54 = reshape(ProjHistMinMax(RUq_min5,'hor',q) + ProjHistMinMax(LDq_min5,'hor',q) + ProjHistMinMax(RDq_min5,'hor',q) + ProjHistMinMax(LUq_min5,'hor',q) + ...
                  ProjHistMinMax(RDq_min5,'ver',q) + ProjHistMinMax(LUq_min5,'ver',q) + ProjHistMinMax(RUq_min5,'ver',q) + ProjHistMinMax(LDq_min5,'ver',q),[],1);
g.max54 = reshape(ProjHistMinMax(RUq_max5,'hor',q) + ProjHistMinMax(LDq_max5,'hor',q) + ProjHistMinMax(RDq_max5,'hor',q) + ProjHistMinMax(LUq_max5,'hor',q) + ...
                  ProjHistMinMax(RDq_max5,'ver',q) + ProjHistMinMax(LUq_max5,'ver',q) + ProjHistMinMax(RUq_max5,'ver',q) + ProjHistMinMax(LDq_max5,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;


end

function g = all3x3(X,q)
% This function outputs co-occurrences of ALL residuals based on the
% KB kernel and its "halves" (EDGE residuals) as listed in Figure 1
% in our journal HUGO paper (version from June 14), including the naming
% convention.
% spam11 (old name KB residual), good, non-directional, hv-symmetrical, to be symmetrized as spam
D = Residual(X,2,'KB');
g.spam11 = reshape(ProjHistSpam(D,'hor',q) + ProjHistSpam(D,'ver',q),[],1); settings.seedIndex=settings.seedIndex+1;
% EDGE residuals
D = Residual(X,2,'edge-h');Du = D(:,1:size(D,2)/2);Db = D(:,size(D,2)/2+1:end);
D = Residual(X,2,'edge-v');Dl = D(:,1:size(D,2)/2);Dr = D(:,size(D,2)/2+1:end);
% spam14h,v  not bad, directional, hv-nonsym, to be symmetrized as spam
g.spam14v = reshape(ProjHistSpam(Du,'ver',q) + ProjHistSpam(Db,'ver',q) + ProjHistSpam(Dl,'hor',q) + ProjHistSpam(Dr,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;
g.spam14h = reshape(ProjHistSpam(Du,'hor',q) + ProjHistSpam(Db,'hor',q) + ProjHistSpam(Dl,'ver',q) + ProjHistSpam(Dr,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax24 -- EXCELLENT, directional, hv-sym, to be symmetrized as mnmx
[Dmin1,Dmin2,Dmin3,Dmin4] = deal(min(Du,Dl),min(Db,Dr),min(Du,Dr),min(Db,Dl));
g.min24 = reshape(ProjHistMinMax(Dmin1,'ver',q) + ProjHistMinMax(Dmin2,'ver',q) + ProjHistMinMax(Dmin3,'ver',q) + ProjHistMinMax(Dmin4,'ver',q) + ... 
                  ProjHistMinMax(Dmin1,'hor',q) + ProjHistMinMax(Dmin2,'hor',q) + ProjHistMinMax(Dmin3,'hor',q) + ProjHistMinMax(Dmin4,'hor',q),[],1);
[Dmax1,Dmax2,Dmax3,Dmax4] = deal(max(Du,Dl),max(Db,Dr),max(Du,Dr),max(Db,Dl));
g.max24 = reshape(ProjHistMinMax(Dmax1,'ver',q) + ProjHistMinMax(Dmax2,'ver',q) + ProjHistMinMax(Dmax3,'ver',q) + ProjHistMinMax(Dmax4,'ver',q) + ... 
                  ProjHistMinMax(Dmax1,'hor',q) + ProjHistMinMax(Dmax2,'hor',q) + ProjHistMinMax(Dmax3,'hor',q) + ProjHistMinMax(Dmax4,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax22 - hv-nonsymmetrical
% min22h -- good, to be symmetrized as mnmx, directional, hv-nonsymmetrical
% min22v -- EXCELLENT - to be symmetrized as mnmx, directional,
[UEq_min,REq_min] = deal(min(Du,Db),min(Dr,Dl));
g.min22h = reshape(ProjHistMinMax(UEq_min,'hor',q) + ProjHistMinMax(REq_min,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
g.min22v = reshape(ProjHistMinMax(UEq_min,'ver',q) + ProjHistMinMax(REq_min,'hor',q),[],1);settings.seedIndex=settings.seedIndex-1;
[UEq_max,REq_max] = deal(max(Du,Db),max(Dr,Dl));
g.max22h = reshape(ProjHistMinMax(UEq_max,'hor',q) + ProjHistMinMax(REq_max,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
g.max22v = reshape(ProjHistMinMax(UEq_max,'ver',q) + ProjHistMinMax(REq_max,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax41 -- good, non-directional, hv-sym, to be symmetrized as mnmx
[Dmin5,Dmax5] = deal(min(Dmin1,Dmin2),max(Dmax1,Dmax2));
g.min41 = reshape(ProjHistMinMax(Dmin5,'ver',q) + ProjHistMinMax(Dmin5,'hor',q),[],1);
g.max41 = reshape(ProjHistMinMax(Dmax5,'ver',q) + ProjHistMinMax(Dmax5,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;
end

function g = all5x5(X,q)
% This function outputs co-occurrences of ALL residuals based on the
% KV kernel and its "halves" (EDGE residuals) as listed in Figure 1
% in our journal HUGO paper (version from June 14), including the naming
% convention.
[M N] = size(X); [I,J] = deal(3:M-2,3:N-2);
% spam11 (old name KV residual), good, non-directional, hv-symmetrical, to be symmetrized as spam
D = Residual(X,3,'KV');
g.spam11 = reshape(ProjHistSpam(D,'hor',q) + ProjHistSpam(D,'ver', q),[],1);settings.seedIndex=settings.seedIndex+1;
% EDGE residuals    
Du = 8*X(I,J-1)+8*X(I-1,J)+8*X(I,J+1)-6*X(I-1,J-1)-6*X(I-1,J+1)-2*X(I,J-2)-2*X(I,J+2)-2*X(I-2,J)+2*X(I-1,J-2)+2*X(I-2,J-1)+2*X(I-2,J+1)+2*X(I-1,J+2)-X(I-2,J-2)-X(I-2,J+2)-12*X(I,J);
Dr = 8*X(I-1,J)+8*X(I,J+1)+8*X(I+1,J)-6*X(I-1,J+1)-6*X(I+1,J+1)-2*X(I-2,J)-2*X(I+2,J)-2*X(I,J+2)+2*X(I-2,J+1)+2*X(I-1,J+2)+2*X(I+1,J+2)+2*X(I+2,J+1)-X(I-2,J+2)-X(I+2,J+2)-12*X(I,J);
Db = 8*X(I,J+1)+8*X(I+1,J)+8*X(I,J-1)-6*X(I+1,J+1)-6*X(I+1,J-1)-2*X(I,J-2)-2*X(I,J+2)-2*X(I+2,J)+2*X(I+1,J+2)+2*X(I+2,J+1)+2*X(I+2,J-1)+2*X(I+1,J-2)-X(I+2,J+2)-X(I+2,J-2)-12*X(I,J);
Dl = 8*X(I+1,J)+8*X(I,J-1)+8*X(I-1,J)-6*X(I+1,J-1)-6*X(I-1,J-1)-2*X(I-2,J)-2*X(I+2,J)-2*X(I,J-2)+2*X(I+2,J-1)+2*X(I+1,J-2)+2*X(I-1,J-2)+2*X(I-2,J-1)-X(I+2,J-2)-X(I-2,J-2)-12*X(I,J);
% spam14v  not bad, directional, hv-nonsym, to be symmetrized as spam
g.spam14v = reshape(ProjHistSpam(Du,'ver', q) + ProjHistSpam(Db,'ver', q) + ProjHistSpam(Dl,'hor',q) + ProjHistSpam(Dr,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;
g.spam14h = reshape(ProjHistSpam(Du,'hor', q) + ProjHistSpam(Db,'hor', q) + ProjHistSpam(Dl,'ver',q) + ProjHistSpam(Dr,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax24 -- EXCELLENT, directional, hv-sym, to be symmetrized as mnmx
[Dmin1,Dmin2,Dmin3,Dmin4] = deal(min(Du,Dl),min(Db,Dr),min(Du,Dr),min(Db,Dl));
g.min24 = reshape(ProjHistMinMax(Dmin1,'ver',q) + ProjHistMinMax(Dmin2,'ver',q) + ProjHistMinMax(Dmin3,'ver',q) + ProjHistMinMax(Dmin4,'ver',q) + ...
                  ProjHistMinMax(Dmin1,'hor',q) + ProjHistMinMax(Dmin2,'hor',q) + ProjHistMinMax(Dmin3,'hor',q) + ProjHistMinMax(Dmin4,'hor',q),[],1);

[Dmax1,Dmax2,Dmax3,Dmax4] = deal(max(Du,Dl),max(Db,Dr),max(Du,Dr),max(Db,Dl));
g.max24 = reshape(ProjHistMinMax(Dmax1,'ver',q) + ProjHistMinMax(Dmax2,'ver',q) + ProjHistMinMax(Dmax3,'ver',q) + ProjHistMinMax(Dmax4,'ver',q) + ...
                  ProjHistMinMax(Dmax1,'hor',q) + ProjHistMinMax(Dmax2,'hor',q) + ProjHistMinMax(Dmax3,'hor',q) + ProjHistMinMax(Dmax4,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax22 - hv-nonsymmetrical
% min22h -- good, to be symmetrized as mnmx, directional, hv-nonsymmetrical
% min22v -- EXCELLENT - to be symmetrized as mnmx, directional,
[UEq_min,REq_min] = deal(min(Du,Db),min(Dr,Dl));
g.min22h = reshape(ProjHistMinMax(UEq_min,'hor',q) + ProjHistMinMax(REq_min,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
g.min22v = reshape(ProjHistMinMax(UEq_min,'ver',q) + ProjHistMinMax(REq_min,'hor',q),[],1);settings.seedIndex=settings.seedIndex-1;
[UEq_max,REq_max] = deal(max(Du,Db),max(Dr,Dl));
g.max22h = reshape(ProjHistMinMax(UEq_max,'hor',q) + ProjHistMinMax(REq_max,'ver',q),[],1);settings.seedIndex=settings.seedIndex+1;
g.max22v = reshape(ProjHistMinMax(UEq_max,'ver',q) + ProjHistMinMax(REq_max,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;
% minmax41 -- good, non-directional, hv-sym, to be symmetrized as mnmx
[Dmin5,Dmax5] = deal(min(Dmin1,Dmin2),max(Dmax1,Dmax2));
g.min41 = reshape(ProjHistMinMax(Dmin5,'ver',q) + ProjHistMinMax(Dmin5,'hor',q),[],1);
g.max41 = reshape(ProjHistMinMax(Dmax5,'ver',q) + ProjHistMinMax(Dmax5,'hor',q),[],1);settings.seedIndex=settings.seedIndex+1;
end

function h = ProjHistSpam(D, type, centerVal)
    RandStream.setGlobalStream(RandStream('mt19937ar','Seed',settings.seedIndex));
    h = zeros(numel(settings.neighborhoods) * settings.qBins * settings.projCount, 1);
    for neighborhoodIndex = 1:numel(settings.neighborhoods)
        neighborhood = settings.neighborhoods{neighborhoodIndex};
        for projIndex = 1:settings.projCount
            P = zeros(size(neighborhood));
            P(neighborhood==1) = randi(settings.B*2+1, sum(neighborhood(:)==1), 1)-settings.B-1;
            if strcmp(type, 'ver'), P = P'; end;
            binSize = sqrt(sum(P(:).^2));
            binEdges = [0:binSize*centerVal:binSize*(settings.qBins-1)*centerVal, inf];
            proj = int32(conv2(D, P, 'valid'));
            h_neigh = histc(abs(proj(:)), binEdges); 
            h_neigh = h_neigh(1:end-1);
            if size(P, 2) > 1
                proj = int32(conv2(D, fliplr(P), 'valid'));
                t = histc(abs(proj(:)), binEdges); 
                h_neigh = h_neigh + t(1:end-1);
            end
            if size(P, 1) > 1
                proj = int32(conv2(D, flipud(P), 'valid'));
                t = histc(abs(proj(:)), binEdges); 
                h_neigh = h_neigh + t(1:end-1);
            end
            if all(size(P)>1)
                proj = int32(conv2(D, rot90(P, 2), 'valid'));
                t = histc(abs(proj(:)), binEdges); 
                h_neigh = h_neigh + t(1:end-1);
            end
            
            h((neighborhoodIndex-1)*settings.qBins*settings.projCount + (projIndex-1)*settings.qBins + 1:(neighborhoodIndex-1)*settings.qBins*settings.projCount + projIndex*settings.qBins, 1) = h_neigh;
        end
    end
end

function h = ProjHistMinMax(D, type, centerVal)
    RandStream.setGlobalStream(RandStream('mt19937ar','Seed',settings.seedIndex));
    h = zeros(numel(settings.neighborhoods) * settings.qBins * settings.projCount * 2, 1);
    for neighborhoodIndex = 1:numel(settings.neighborhoods)
        neighborhood = settings.neighborhoods{neighborhoodIndex};
        for projIndex = 1:settings.projCount
            P = zeros(size(neighborhood));
            P(neighborhood==1) = randi(settings.B*2+1, sum(neighborhood(:)==1), 1)-settings.B-1;
            if strcmp(type, 'ver'), P = P'; end;
            binSize = sqrt(sum(P(:).^2));
            binEdges = [-inf, -binSize*(settings.qBins-1)*centerVal:binSize*centerVal:binSize*(settings.qBins-1)*centerVal, inf];
            proj = int32(conv2(D, P, 'valid'));
            h_neigh = histc(proj(:), binEdges); 
            h_neigh = h_neigh(1:end-1);
            if size(P, 2) > 1
                proj = int32(conv2(D, fliplr(P), 'valid'));
                t = histc(proj(:), binEdges); 
                h_neigh = h_neigh + t(1:end-1);
            end
            if size(P, 1) > 1
                proj = int32(conv2(D, flipud(P), 'valid'));
                t = histc(proj(:), binEdges); 
                h_neigh = h_neigh + t(1:end-1);
            end
            if all(size(P)>1)
                proj = int32(conv2(D, rot90(P, 2), 'valid'));
                t = histc(proj(:), binEdges); 
                h_neigh = h_neigh + t(1:end-1);
            end
            
            h((neighborhoodIndex-1)*settings.qBins*2*settings.projCount + (projIndex-1)*settings.qBins*2 + 1:(neighborhoodIndex-1)*settings.qBins*2*settings.projCount + projIndex*settings.qBins*2, 1) = h_neigh;
        end
    end
 end

function result = mergeMinMax(Fmin, Fmax)
    Fmin = reshape(Fmin, 2*settings.qBins, []);
    Fmin = flipud(Fmin);
    result = Fmax + Fmin(:)';
end

function D = Residual(X,order,type)
% Computes the noise residual of a given type and order from MxN image X.
% residual order \in {1,2,3,4,5,6}
% type \in {hor,ver,diag,mdiag,KB,edge-h,edge-v,edge-d,edge-m}
% The resulting residual is an (M-b)x(N-b) array of the specified order,
% where b = ceil(order/2). This cropping is little more than it needs to 
% be to make sure all the residuals are easily "synchronized".
% !!!!!!!!!!!!! Use order = 2 with KB and all edge residuals !!!!!!!!!!!!!

[M N] = size(X);
I = 1+ceil(order/2) : M-ceil(order/2);
J = 1+ceil(order/2) : N-ceil(order/2);

switch type
    case 'hor'
        switch order
            case 1, D = - X(I,J) + X(I,J+1);
            case 2, D = X(I,J-1) - 2*X(I,J) + X(I,J+1);
            case 3, D = X(I,J-1) - 3*X(I,J) + 3*X(I,J+1) - X(I,J+2);
            case 4, D = -X(I,J-2) + 4*X(I,J-1) - 6*X(I,J) + 4*X(I,J+1) - X(I,J+2);
            case 5, D = -X(I,J-2) + 5*X(I,J-1) - 10*X(I,J) + 10*X(I,J+1) - 5*X(I,J+2) + X(I,J+3);
            case 6, D = X(I,J-3) - 6*X(I,J-2) + 15*X(I,J-1) - 20*X(I,J) + 15*X(I,J+1) - 6*X(I,J+2) + X(I,J+3);
        end
    case 'ver'
        switch order
            case 1, D = - X(I,J) + X(I+1,J);
            case 2, D = X(I-1,J) - 2*X(I,J) + X(I+1,J);
            case 3, D = X(I-1,J) - 3*X(I,J) + 3*X(I+1,J) - X(I+2,J);
            case 4, D = -X(I-2,J) + 4*X(I-1,J) - 6*X(I,J) + 4*X(I+1,J) - X(I+2,J);
            case 5, D = -X(I-2,J) + 5*X(I-1,J) - 10*X(I,J) + 10*X(I+1,J) - 5*X(I+2,J) + X(I+3,J);
            case 6, D = X(I-3,J) - 6*X(I-2,J) + 15*X(I-1,J) - 20*X(I,J) + 15*X(I+1,J) - 6*X(I+2,J) + X(I+3,J);
        end
    case 'diag'
        switch order
            case 1, D = - X(I,J) + X(I+1,J+1);
            case 2, D = X(I-1,J-1) - 2*X(I,J) + X(I+1,J+1);
            case 3, D = X(I-1,J-1) - 3*X(I,J) + 3*X(I+1,J+1) - X(I+2,J+2);
            case 4, D = -X(I-2,J-2) + 4*X(I-1,J-1) - 6*X(I,J) + 4*X(I+1,J+1) - X(I+2,J+2);
            case 5, D = -X(I-2,J-2) + 5*X(I-1,J-1) - 10*X(I,J) + 10*X(I+1,J+1) - 5*X(I+2,J+2) + X(I+3,J+3);
            case 6, D = X(I-3,J-3) - 6*X(I-2,J-2) + 15*X(I-1,J-1) - 20*X(I,J) + 15*X(I+1,J+1) - 6*X(I+2,J+2) + X(I+3,J+3);
        end
    case 'mdiag'
        switch order
            case 1, D = - X(I,J) + X(I-1,J+1);
            case 2, D = X(I-1,J+1) - 2*X(I,J) + X(I+1,J-1);
            case 3, D = X(I-1,J+1) - 3*X(I,J) + 3*X(I+1,J-1) - X(I+2,J-2);
            case 4, D = -X(I-2,J+2) + 4*X(I-1,J+1) - 6*X(I,J) + 4*X(I+1,J-1) - X(I+2,J-2);
            case 5, D = -X(I-2,J+2) + 5*X(I-1,J+1) - 10*X(I,J) + 10*X(I+1,J-1) - 5*X(I+2,J-2) + X(I+3,J-3);
            case 6, D = X(I-3,J+3) - 6*X(I-2,J+2) + 15*X(I-1,J+1) - 20*X(I,J) + 15*X(I+1,J-1) - 6*X(I+2,J-2) + X(I+3,J-3);
        end
    case 'KB'
        D = -X(I-1,J-1) + 2*X(I-1,J) - X(I-1,J+1) + 2*X(I,J-1) - 4*X(I,J) + 2*X(I,J+1) - X(I+1,J-1) + 2*X(I+1,J) - X(I+1,J+1);
    case 'edge-h'
        Du = 2*X(I-1,J) + 2*X(I,J-1) + 2*X(I,J+1) - X(I-1,J-1) - X(I-1,J+1) - 4*X(I,J);   %   -1  2 -1
        Db = 2*X(I+1,J) + 2*X(I,J-1) + 2*X(I,J+1) - X(I+1,J-1) - X(I+1,J+1) - 4*X(I,J);   %    2  C  2    +  flipped vertically
        D = [Du,Db];
    case 'edge-v'
        Dl = 2*X(I,J-1) + 2*X(I-1,J) + 2*X(I+1,J) - X(I-1,J-1) - X(I+1,J-1) - 4*X(I,J);   %   -1  2
        Dr = 2*X(I,J+1) + 2*X(I-1,J) + 2*X(I+1,J) - X(I-1,J+1) - X(I+1,J+1) - 4*X(I,J);   %    2  C       +  flipped horizontally
        D = [Dl,Dr];                                                                      %   -1  2
    case 'edge-m'
        Dlu = 2*X(I,J-1) + 2*X(I-1,J) - X(I-1,J-1) - X(I+1,J-1) - X(I-1,J+1) - X(I,J); %      -1  2 -1
        Drb = 2*X(I,J+1) + 2*X(I+1,J) - X(I+1,J+1) - X(I+1,J-1) - X(I-1,J+1) - X(I,J); %       2  C       +  flipped mdiag
        D = [Dlu,Drb];                                                                 %      -1
    case 'edge-d'
        Dru = 2*X(I-1,J) + 2*X(I,J+1) - X(I-1,J+1) - X(I-1,J-1) - X(I+1,J+1) - X(I,J); %      -1  2 -1
        Dlb = 2*X(I,J-1) + 2*X(I+1,J) - X(I+1,J-1) - X(I+1,J+1) - X(I-1,J-1) - X(I,J); %          C  2    +  flipped diag
        D = [Dru,Dlb];                                                                 %            -1
    case 'KV'
        D = 8*X(I-1,J) + 8*X(I+1,J) + 8*X(I,J-1) + 8*X(I,J+1);
        D = D - 6*X(I-1,J+1) - 6*X(I-1,J-1) - 6*X(I+1,J-1) - 6*X(I+1,J+1);
        D = D - 2*X(I-2,J) - 2*X(I+2,J) - 2*X(I,J+2) - 2*X(I,J-2);
        D = D + 2*X(I-1,J-2) + 2*X(I-2,J-1) + 2*X(I-2,J+1) + 2*X(I-1,J+2) + 2*X(I+1,J+2) + 2*X(I+2,J+1) + 2*X(I+2,J-1) + 2*X(I+1,J-2);
        D = D - X(I-2,J-2) - X(I-2,J+2) - X(I+2,J-2) - X(I+2,J+2) - 12*X(I,J);
end

end

end