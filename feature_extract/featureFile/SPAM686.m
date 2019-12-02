function F = SPAM686(IMG)
% -------------------------------------------------------------------------
% Contact: jan@kodovsky.com | fridrich@binghamton.edu | November 2011
%          http://dde.binghamton.edu/download/feature_extractors
% -------------------------------------------------------------------------
% Extracts spatial domain SPAM features, 2nd order, T=3. Dimensionality
% 686. See the original publication [1] for more details. In this
% implementation, we incorporated the following two modifications as
% compared to [1] (slightly better performance):
%  - marginalization of elements out of [-T,T] into borders (NOT throwing
%    them away)
%  - joint sample distribution (co-occurrence) instead of transition
%    probability matrix (different normalization)
% -------------------------------------------------------------------------
% Input: IMG ... input images (can be also JPEG image)
% Output: F .... resulting SPAM features
% -------------------------------------------------------------------------
% [1] T. Pevny, P. Bas, and J. Fridrich, Steganalysis by Subtractive Pixel
% Adjacency Matrix IEEE Trans. on Info. Forensics and Security, vol. 5(2),
% pp. 215鈥?24, 2010.
% -------------------------------------------------------------------------

% 若IMG是像素矩阵,则调用
if(isnumeric(IMG))
    F = spam_extract_2( double(IMG),3 );
else
    F = spam_extract_2(double(imread(IMG)),3);
end
F = F';


function F = spam_extract_2(X,T)

% horizontal left-right
D = X(:,1:end-1) - X(:,2:end);
L = D(:,3:end); C = D(:,2:end-1); R = D(:,1:end-2);
Mh1 = GetM3(L,C,R,T);

% horizontal right-left
D = -D;
L = D(:,1:end-2); C = D(:,2:end-1); R = D(:,3:end);
Mh2 = GetM3(L,C,R,T);
% vertical bottom top
D = X(1:end-1,:) - X(2:end,:);
L = D(3:end,:); C = D(2:end-1,:); R = D(1:end-2,:);
Mv1 = GetM3(L,C,R,T);

% vertical top bottom
D = -D;
L = D(1:end-2,:); C = D(2:end-1,:); R = D(3:end,:);
Mv2 = GetM3(L,C,R,T);

% diagonal left-right
D = X(1:end-1,1:end-1) - X(2:end,2:end);
L = D(3:end,3:end); C = D(2:end-1,2:end-1); R = D(1:end-2,1:end-2);
Md1 = GetM3(L,C,R,T);

% diagonal right-left
D = -D;
L = D(1:end-2,1:end-2); C = D(2:end-1,2:end-1); R = D(3:end,3:end);
Md2 = GetM3(L,C,R,T);

% minor diagonal left-right
D = X(2:end,1:end-1) - X(1:end-1,2:end);
L = D(1:end-2,3:end); C = D(2:end-1,2:end-1); R = D(3:end,1:end-2);
Mm1 = GetM3(L,C,R,T);

% minor diagonal right-left
D = -D;
L = D(3:end,1:end-2); C = D(2:end-1,2:end-1); R = D(1:end-2,3:end);
Mm2 = GetM3(L,C,R,T);

F1 = (Mh1+Mh2+Mv1+Mv2)/4;
F2 = (Md1+Md2+Mm1+Mm2)/4;
F = [F1;F2];

function M = GetM3(L,C,R,T)
% marginalization into borders
L = L(:); L(L<-T) = -T; L(L>T) = T;
C = C(:); C(C<-T) = -T; C(C>T) = T;
R = R(:); R(R<-T) = -T; R(R>T) = T;

% get cooccurences [-T...T]
M = zeros(2*T+1,2*T+1,2*T+1);
for i=-T:T
    C2 = C(L==i);
    R2 = R(L==i);
    for j=-T:T
        R3 = R2(C2==j);
        for k=-T:T
            M(i+T+1,j+T+1,k+T+1) = sum(R3==k);
        end
    end
end

% normalization
M = M(:)/sum(M(:));