%% Find the complex in a multivariate autoregressive model, 
%   X^t = A*X^{t-1} + E,
% where A is a connectivity matrix and E is gaussian noise.

%% options
%           options: options for computing phi, MIP search, and complex
%           search
%           
%           options.type_of_dist:
%              'Gauss': Gaussian distribution
%              'discrete': discrete probability distribution
%              'UndirectedGraph': Undirected Graph
%           options.type_of_phi:
%              'MI1': Multi (Mutual) information. I(X_1; X_2). (IIT1.0)
%              'MI': Multi (Mutual) information, I(X_1, Y_1; X_2, Y_2)
%              'SI': phi_H, stochastic interaction
%              'star': phi_star, based on mismatched decoding
%              'Geo': phi_G, information geometry version
%           options.type_of_MIPsearch
%              'Exhaustive': exhaustive search
%              'Queyranne': Queyranne algorithm
%              'REMCMC': Replica Exchange Monte Carlo Method 
%              'StoerWagner': mincut search algorithm for undirected graphs
%           options.type_of_complexsearch
%               'Exhaustive': exhaustive search
%               'Recursive': recursive MIP search
%           options.normalization
%              (Available only when the dist. is 'Gauss' and the complex search is 'Exhaustive')
%              Regarding normalization of phi by Entropy when searching
%              the MIPs. 
%                 0: without normalization (default)
%                 1: with normalization
%              Note that, after finding the MIPs, phi wo/ normalization at
%              the MIPs is used to compare subsets and find the complexes in both options. 

clear all;
addpath(genpath('../PhiToolbox_Graph'))


%% generate data
disp('Generating data...')

%%% construct connectivity matrix %%%
N = 7; % the number of elements
A = eye(N)/N; % A: connectivity matrix
A(3:N-2, 3:N-2) = 1/N;
A(1:2, 1:2) = 1/N;
A(6:7, 6:7) = 0.5/N;
A = A + 0.01*randn(N)/sqrt(N);

%%% generate time series X using the AR model %%%
T = 10^6;
X = zeros(N,T);
X(:,1) = randn(N,1);
for t=2: T
    E = randn(N,1);
    X(:,t) = A*X(:,t-1) + E;
end


%% Convert data to an undirected graph using pairwise phi
params_ToGenerateGraph.tau = 1;

options_ToGenerateGraph.type_of_dist = 'Gauss';
options_ToGenerateGraph.type_of_phi = 'MI1';

%%% Compute pairwise phi and set it as the weight of the graph %%%
g = zeros(N,N);
Z = [1 2];
for i=1:N
    for j = (i+1):N     
        g(i,j) = phi_comp(X([i j],:),Z,params_ToGenerateGraph,options_ToGenerateGraph);
        g(j,i) = g(i,j);
    end
end
G = graph(g);

figure(1)
subplot(1,2,1), imagesc(A)
title('Connectivity matrix of AR model')

subplot(1,2,2), imagesc(g)
title('Undirected graph')

drawnow


%% params
params.tau = 1; % time lag

%% options
probs.g = g;
probs.number_of_elements = size(probs.g,1);

options.type_of_dist = 'UndirectedGraph';
options.type_of_MIPsearch = 'StoerWagner';

tic;
%% Recursive
[complexes, phis_complexes, main_complexes, phis_main_complexes, Res] = ...
    Complex_Recursive( probs, options );

t = toc;

%% Show results
main_complexes_str = cell(size(main_complexes));
for i = 1:length(main_complexes)
    main_complexes_str{i} =  num2str(main_complexes{i});
end
figure
bar(phis_main_complexes)
set(gca, 'xticklabel', main_complexes_str)
title('Main Complexes')
ylabel('$w^\mathrm{mc}$', 'Interpreter', 'latex', 'FontSize', 18)
xlabel('Indices of the main complexes')


thetas = 2*pi*(0:(N-1))/N;
XCoor = cos(thetas);
YCoor = sin(thetas);

figure
colormap(KovesiRainbow)
plotComplexes(complexes,phis_complexes,'BirdsEye','Graph', G, 'XData',XCoor,'YData',YCoor, 'LineWidth', 2)
title('(Main) Complexes')
colorbar

min_phi_comp = min(phis_complexes);
max_phi_comp = max(phis_complexes);
zlim([min_phi_comp, max_phi_comp])
zlabel('$w^\mathrm{mc}$', 'Interpreter', 'latex', 'FontSize', 18)