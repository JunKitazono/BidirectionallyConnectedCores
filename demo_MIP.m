
% Find Minimum Information Partition (MIP) in a multivariate autoregressive model, 
%   X^t = A*X^{t-1} + E,
% where A is a connectivity matrix and E is gaussian noise.

%% generate data
disp('Generating data...')

%%% construct connectivity matrix %%%
N = 10; % the number of elements
A = zeros(N); % A: connectivity matrix (block structured)
n_block = 2;
for i = 1:n_block
    indices_block = ((1+(i-1)*N/n_block):i*N/n_block)';
    A(indices_block, indices_block) = 1/N;
end
A = A + 0.01*randn(N)/sqrt(N);

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

%%% estimate covariance matrices %%%
% Note: If the sample size is small and the number of elements is large, the
% empirical estimators of covariance matrices below will be very unstable.
% see https://en.wikipedia.org/wiki/Estimation_of_covariance_matrices.

tau = 1; % time lag
t_range1 = 1: 1: T-tau;
t_range2 = 1+tau: 1: T;

X1 = X(:,t_range1);
X1 = bsxfun(@minus, X1, mean(X1,2)); % subtract mean
X2 = X(:,t_range2);
X2 = bsxfun(@minus, X2, mean(X2,2)); % subtract mean

Cov_X = X1*X1'/(T-tau-1); % equal-time covariance matrix at "PAST"
Cov_Y = X2*X2'/(T-tau-1); % equal-time covariance matrix at "PRESENT"
Cov_XY = X1*X2'/(T-tau-1); % time-lagged covariance matrix


%% find Minimum Information Partition (MIP)

% type_of_phi:
%    'MI1': Multi (Mutual) information, e.g., I(X_1; X_2). (IIT1.0)
%    'MI': Multi (Mutual) information, e.g., I(X_1, Y_1; X_2, Y_2)
%    'SI': phi_H, stochastic interaction
%    'star': phi_star, based on mismatched decoding
%    'Geo': phi_G, information geometry version
type_of_phi = 'SI';
type_of_dist = 'Gauss';

%%% Exhaustive Search %%%
disp('Exhaustive Search...')

%% pre-compute covariances
tic;
probs{1} = Cov_X;
probs{2} = Cov_XY;
probs{3} = Cov_Y;
params(1) = N;
params(2) = tau;

[Z_MIP, phi_MIP] = MIP_Exhaustive_probs( type_of_dist, type_of_phi, params, probs);

t_Exhaustive = toc;
disp('Exhaustive Search finished.')
disp(['phi at the MIP: ', num2str(phi_MIP)])
disp(['the MIP: ', num2str(Z_MIP)])

%% without covariances
tic;

[Z_MIP, phi_MIP] = MIP_Exhaustive( type_of_dist, type_of_phi, X, params);
t_Exhaustive = toc;
disp('Exhaustive Search finished.')
disp(['phi at the MIP: ', num2str(phi_MIP)])
disp(['the MIP: ', num2str(Z_MIP)])

%%% Queyeranne's algorithm %%%
% disp('Queyranne''s algorithm...')
% tic;
% [Z_MIP_Q, phi_MIP_Q] = MIP_Queyranne( type_of_phi, Cov_X, Cov_XY, Cov_Y );
% t_Queyeranne = toc;
% disp('Queyranne''s algorithm finished.')
% disp(['phi at the estimated MIP: ', num2str(phi_MIP_Q)])
% disp(['the estimated MIP: ', num2str(Z_MIP_Q)])
% 
% %%% Replica Exchange Markov Chain Monte Carlo (REMCMC) %%%
% options = [];
% disp('REMCMC...')
% tic;
% [Z_MIP_REMCMC, phi_MIP_REMCMC, ...
%     Energy_history, State_history, Exchange_history, T_history, wasConverged, NumCalls] = ...
%     MIP_REMCMC( type_of_phi, options, Cov_X, Cov_XY, Cov_Y );
% t_REMCMC = toc;
% disp('REMCMC finished.')
% disp(['phi at the estimated MIP: ', num2str(phi_MIP_REMCMC)])
% disp(['the estimated MIP: ', num2str(Z_MIP_REMCMC)])
% 
