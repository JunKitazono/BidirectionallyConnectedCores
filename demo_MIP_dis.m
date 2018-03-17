addpath(genpath('../PhiToolbox'))

N = 4; % number of units
T = 10^6; % number of iterations
tau = 1; % time delay
N_st = 2;  % number of states

%% generate time series from Boltzman machine
beta = 4; % inverse temperature

% connectivity matrix
W = zeros(N,N);
community = [1 2 2 1];
for i=1: N
    for j=1: N
        if i~=j
            if community(i) == community(j)
                % W(i,j) = 0.2; % for N = 8
                W(i,j) = 0.4;
            else
                W(i,j) = 0;
            end
        else
            W(i,i) = 0.1;
        end
    end
end
X = generate_Boltzmann(beta,W,N,T); 

%% 

T_seg = 1000;
figure(1)
t_vec1 = 1: T_seg;
t_vec2 = 2*10^3: 2*10^3+T_seg;
t_vec3 = 10^4: 10^4+T_seg;
t_vec4 = 10^5: 10^5+T_seg;
t_vec5 = T-300: T;
subplot(3,2,1),imagesc(X(:,t_vec1));
subplot(3,2,2),imagesc(X(:,t_vec2));
subplot(3,2,3),imagesc(X(:,t_vec3));
subplot(3,2,4),imagesc(X(:,t_vec4));
subplot(3,2,5),imagesc(X(:,t_vec5));

%% compute correlation 
R = corrcoef(X');
disp('Correlation Matrix')
disp(R);

%% find Minimum Information Partition (MIP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
type_of_dist = 'discrete';
type_of_phi = 'MI1';

%%%%%%%%%% with pre-computed probability distributions %%%%%%%%%% 
disp('Find the MIP with pre-computed probability distributions')

%% estimate probability distributions
disp('Estimating probability distributions...')
probs = data_to_probs(type_of_dist, X, tau, N_st);

%% Exhaustive Search %%
disp('Exhaustive Search...')
tic;
[Z_MIP_with, phi_MIP_with] = MIP_Exhaustive_probs( type_of_phi, probs );
t_Exhaustive_with = toc;
disp( ['Exhaustive Search finished. CalcTime=', num2str(t_Exhaustive_with)])
disp(['phi at the MIP: ', num2str(phi_MIP_with)])
disp(['the MIP: ', num2str(Z_MIP_with)])
disp(' ')

%% Queyeranne's algorithm %%%
disp('Queyranne''s algorithm...')
tic;
[Z_MIP_Q_with, phi_MIP_Q_with] = MIP_Queyranne_probs( type_of_phi, probs );
t_Queyranne_with = toc;
disp(['Queyranne''s algorithm finished. CalcTime=', num2str(t_Queyranne_with)])
disp(['phi at the estimated MIP: ', num2str(phi_MIP_Q_with)])
disp(['the estimated MIP: ', num2str(Z_MIP_Q_with)])
disp(' ')

%% Replica Exchange Markov Chain Monte Carlo (REMCMC) %%
options.ShowFig = 0;
options.nMCS = 100;
disp('REMCMC...')
tic;
[Z_MIP_REMCMC_with, phi_MIP_REMCMC_with, ...
    phi_history, State_history, Exchange_history, T_history, wasConverged, NumCalls] = ...
    MIP_REMCMC_probs( type_of_phi, probs, options );
t_REMCMC_with = toc;
disp(['REMCMC finished. CalcTime=', num2str(t_REMCMC_with)])
disp(['phi at the estimated MIP: ', num2str(phi_MIP_REMCMC_with)])
disp(['the estimated MIP: ', num2str(Z_MIP_REMCMC_with)])
disp(' ')


%%%%%%%%%% without pre-computed probability distributions %%%%%%%%%%
disp('Find the MIP without pre-computed probability distributions')
%% Exhaustive search %%
disp('Exhaustive Search...')
tic;
[Z_MIP_without, phi_MIP_without] = MIP_Exhaustive( type_of_dist, type_of_phi, X, tau, N_st );
t_Exhaustive_without = toc;
disp(['Exhaustive Search finished. CalcTime=', num2str(t_Exhaustive_without)] )
disp(['phi at the MIP: ', num2str(phi_MIP_without)])
disp(['the MIP: ', num2str(Z_MIP_without)])
disp(' ')

%% Queyranne's algorithm %%
disp('Queyranne''s algorithm...')
tic;
[Z_MIP_Q_without, phi_MIP_Q_without] = MIP_Queyranne( type_of_dist, type_of_phi, X, tau, N_st );
t_Queyranne_without = toc;
disp(['Exhaustive Search finished. CalcTime=', num2str(t_Queyranne_without)])
disp(['phi at the MIP: ', num2str(phi_MIP_Q_without)])
disp(['the MIP: ', num2str(Z_MIP_Q_without)])
disp(' ')

%% Replica Exchange Markov Chain Monte Carlo (REMCMC) %%
options.ShowFig = 0;
options.nMCS = 100;
disp('REMCMC...')
tic;
[Z_MIP_REMCMC_without, phi_MIP_REMCMC_without, ...
    phi_history, State_history, Exchange_history, T_history, wasConverged, NumCalls] = ...
    MIP_REMCMC( type_of_dist, type_of_phi, X, tau, N_st, options );
t_REMCMC_without = toc;
disp(['REMCMC finished. CalcTime=', num2str(t_REMCMC_without)])
disp(['phi at the estimated MIP: ', num2str(phi_MIP_REMCMC_without)])
disp(['the estimated MIP: ', num2str(Z_MIP_REMCMC_without)])
disp(' ')


