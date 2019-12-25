function rho = HILL(cover)
F = [-0.25, 0.5, -0.25; 0.5, -1, 0.5; -0.25, 0.5, -0.25];
%% Get embedding costs
% inicialization
cover = double( cover );

wetCost = 10^10;

% compute residual
R = imfilter(cover, F, 'symmetric', 'conv', 'same');

% compute suitability
xi = imfilter(abs(R), (1/9)*ones(3, 3), 'symmetric', 'conv', 'same');

% compute embedding costs \rho
rho = imfilter(1./xi, (1/225)*ones(15, 15), 'symmetric', 'conv', 'same');
% figure;
% imshow(rho);

% adjust embedding costs
rho(rho > wetCost) = wetCost; % threshold on the costs
rho(isnan(rho)) = wetCost; % if all xi{} are zero threshold the cost

end
