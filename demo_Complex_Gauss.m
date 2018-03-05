
% Find Minimum Information Partition (MIP) in a multivariate autoregressive model, 
%   X^t = A*X^{t-1} + E,
% where A is a connectivity matrix and E is gaussian noise.


%% generate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Generating data...')

%%% construct connectivity matrix %%%
N = 6; % the number of elements
A = zeros(N); % A: connectivity matrix (block structured)
n_block = 2;
A = eye(N)/N;
A(2:N-2, 2:N-2) = 1/N;

h_A = figure;
imagesc(A)
drawnow
title('Connectivity Matrix A')

%%% generate random gaussian time series X %%%
T = 10^6;
X = zeros(N,T);
X(:,1) = randn(N,1);
for t=2: T
    E = randn(N,1);
    X(:,t) = A*X(:,t-1) + E;
end


%% find the complex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type_of_dist = 'Gauss';
% type_of_dist:
%    'Gauss': Gaussian distribution
%    'discrete': discrete probability distribution

type_of_phi = 'MI1';
% type_of_phi:
%    'MI1': Multi (Mutual) information, e.g., I(X_1; X_2). (IIT1.0)
%    'MI': Multi (Mutual) information, e.g., I(X_1, Y_1; X_2, Y_2)
%    'SI': phi_H, stochastic interaction
%    'star': phi_star, based on mismatched decoding
%    'Geo': phi_G, information geometry version

type_of_MIPsearch = 'Queyranne';
% type_of_MIPsearch: 
%    'exhaustive': 
%    'Queyranne': 
%    'REMCMC':

tau = 1; % time lag

% Convert data to covariance matrices
% probs = data_to_probs(type_of_dist, X, tau);

%%% with pre-computed covariances %%%
% [indices_Complex, phi_Complex, indices, phis, Zs] = ...
%     Complex_Exhaustive_probs( type_of_phi, type_of_MIPsearch, probs );

% %%% without pre-computed covariances %%%
[ indices_Complex, phi_Complex, indices, phis, Zs ] = ...
   Complex_Exhaustive( type_of_dist, type_of_phi, type_of_MIPsearch, X, tau);

