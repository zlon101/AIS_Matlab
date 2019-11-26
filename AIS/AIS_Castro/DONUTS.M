function S = donuts;

%
% Ph.D. Thesis
% Leandro Nunes de Castro
% February, 2000
% Function generates the dataset for the 2-DONUTS problem
%

% DONUT 1
step = 0.2;
t = 0:step:4*pi;
l = length(t);
h = 0; r = 1;
x = r .* cos(t);
y = r .* sin(t);
z = h .* ones(1,l);
m1 = [x',y',z'];
% size(x), size(y), size(z), pause
plot3(x',y',z','b*'); hold on;
xlabel('x'); ylabel('y'); zlabel('z');
% axis([-1.2 1.2 -1.2 1.2 -1.2 1.2]);

% DONUT 2
step = 0.2;
t = 0:step:4*pi;
l = length(t);
h = 1; r = 1;
x = h .* ones(1,l);
y = r .* cos(t);
z = r .* sin(t);
m2 = [x',y',z'];
% if m1-m2 < 0.1, disp('menor');end;
% size(x), size(y), size(z), pause
% figure(2); 
plot3(x',y',z','ro'); hold off;
S = [m1;m2];
disp(sprintf('Dataset size: [%d,%d]',size(S)));
