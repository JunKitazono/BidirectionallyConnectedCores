% compute several measures of integrated information in a multivariate autoregressive model, 
% X^t = A*X^{t-1} + E,
% where A is a connectivity matrix and E is gaussian noise.

addpath(genpath('../PhiToolbox'))

N = 2; % the number of elements
A = [0.2 0.1; 0.5 0.2]; % connectivity matrix
Cov_E = eye(N,N); % covariance matrix of E

%% generate random gaussian time series X
T = 10^5;
X = zeros(N,T);
X(:,1) = randn(N,1);
for t=2: T
    E = randn(N,1);
    X(:,t) = A*X(:,t-1) + E;
end

%% compute phi from time series 
type_of_dist = 'Gauss';
% type_of_dist:
%    'Gauss': Gaussian distribution
%    'discrete': discrete probability distribution

% available options of type_of_phi for Gaussian distributions:
%    'MI1': Multi (Mutual) information, e.g., I(X_1; X_2). (IIT1.0)
%    'SI': phi_H, stochastic interaction
%    'star': phi_star, based on mismatched decoding
%    'Geo': phi_G, information geometry version

Z = [1 2]; % partition with which phi is computed

% mutual information
type_of_phi = 'MI1';
MI = phi_comp(type_of_dist, type_of_phi, Z, X, tau);

% stochastic interaction
type_of_phi = 'SI';
SI = phi_comp(type_of_dist, type_of_phi, Z, X, tau);

% phi*
type_of_phi = 'star';
phi_star = phi_comp(type_of_dist, type_of_phi, Z, X, tau);

% geometric phi
type_of_phi = 'Geo';
phi_G = phi_comp(type_of_dist, type_of_phi, Z, X, tau);

%%
fprintf('MI=%f SI=%f phi*=%f phi_G=%f\n',MI,SI,phi_star,phi_G);

%% compute phi from pre-computed covariance matrices

%% estimate covariance matrices 
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

Cov_X = X1*X1'/(T-tau-1); % equal-time covariance matrix at past
Cov_Y = X2*X2'/(T-tau-1); % equal-time covariance matrix at present
Cov_XY = X1*X2'/(T-tau-1); % time-lagged covariance matrix

probs.Cov_X = Cov_X;
probs.Cov_Y = Cov_Y;
probs.Cov_XY = Cov_XY;

% the above method of estimating covariance matrices from time series X
% is implemented in the function "data_to_probs"
% probs = data_to_probs(type_of_dist, X, tau); 

%% compute phi
Z = [1 2]; % labels of subsystems

% mutual information
type_of_phi = 'MI1';
MI = phi_comp_probs(type_of_dist, type_of_phi, Z, probs);

% stochastic interaction
type_of_phi = 'SI';
SI = phi_comp_probs(type_of_dist, type_of_phi, Z, probs);

% phi*
type_of_phi = 'star';
phi_star = phi_comp_probs(type_of_dist, type_of_phi, Z, probs);

% geometric phi
type_of_phi = 'Geo';
phi_G = phi_comp_probs(type_of_dist, type_of_phi, Z, probs);

%%
fprintf('MI=%f SI=%f phi*=%f phi_G=%f\n',MI,SI,phi_star,phi_G);