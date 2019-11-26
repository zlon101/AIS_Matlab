function S = spir;

%
% Ph.D. Thesis
% Leandro Nunes de Castro
% February, 2000
% Function generates the dataset for the SPIRAL problem
%

% SPIRAL 1
step = 0.1;
t = 0:step:3*pi;
l = length(t);
li = 0.1; ls = 1;
r = li:(ls-li)/(l-1):ls;
% aux = 0.02 .* randn(1,l);
aux = 1;
x = aux + (r+.1) .* cos(t);
y = aux + (r+.1) .* sin(t);
m1 = [x',y'];
% size(m1),
plot(x',y','b*');

% SPIRAL 2
t = 0:step:3*pi;
l = length(t);
li = 0; ls = 1;
r = li:(ls-li)/(l-1):ls;
% aux = 0.02 .* randn(1,l);
aux = 1;
x = aux - (r+.2) .* cos(t);
y = aux - (r+.2) .* sin(t);
m2 = [x',y'];
hold on; plot(x',y','ro');
hold off;
S = [m1;m2];
disp(sprintf('Dataset size: [%d,%d]',size(S)));
