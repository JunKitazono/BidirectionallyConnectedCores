
% Find the complex.
% Monkey ECoG data (Neurotycho.org)
% '20120730PF_Anesthesia+and+Sleep_Chibi_Toru+Yanagawa_mat_ECoG128'
addpath(genpath('../PhiToolbox'))

%% load datasets

condition = 'awake'; 
% condition = 'anes';

% extract 1-minute signal
window_length = 60*1000; % 1 minute
subsampling_freq = 10; % Down-sample from 1kHz to 100Hz

switch condition
    case 'awake'
        load('Neurotycho/Data_AwakeEyesClosed.mat')
        X = X_AwakeEyesClosed(:, 1:subsampling_freq:window_length);
        % X_AwakeEyesClosed: 9 minutes signals of 64 channeles. 64 by 54000 (=9 minutes * 60 sec. * 1000Hz) matrix.
    case 'anes'
        load('Neurotycho/Data_Anesthetized.mat')
        X = X_Anesthetized(:, 1:subsampling_freq:window_length);
        % X_Anesthetized: 9 minutes signals of 64 channeles. 64 by 54000 (=9 minutes* 60 sec. * 1000Hz) matrix.
end

%% 
target_ch = 1:62;
N = length(target_ch);
X = X(target_ch,:);


%% find the complex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.type_of_dist = 'Gauss';
% type_of_dist:
%    'Gauss': Gaussian distribution
%    'discrete': discrete probability distribution

options.type_of_phi = 'MI1';
% type_of_phi:
%    'MI1': Multi (Mutual) information, e.g., I(X_1; X_2). (IIT1.0)
%    'MI': Multi (Mutual) information, e.g., I(X_1, Y_1; X_2, Y_2)
%    'SI': phi_H, stochastic interaction
%    'star': phi_star, based on mismatched decoding
%    'Geo': phi_G, information geometry version

options.type_of_MIPsearch = 'Queyranne';
% type_of_MIPsearch: 
%    'exhaustive': 
%    'Queyranne': 
%    'REMCMC':

params.tau = 1; % time lag

probs = data_to_probs( X, params, options );

[complexes, phis_complex, Res] = Complex_Recursive_probs( probs, options );

[phi_largest, row_phi_largest] = max(phis_complex);
complex_phi_largest = complexes{row_phi_largest};

%% make Figures
CortexMap = load('ChibiMap_bipolar.mat');

% plot the comlex with largest phi
figure(1)
imagesc(CortexMap.I*(1/3)+2*256/3), axis equal
hold on
scatter(CortexMap.X(target_ch), CortexMap.Y(target_ch), 'r')
scatter(CortexMap.X(complex_phi_largest), CortexMap.Y(complex_phi_largest), 'r', 'filled') 
title('Complex with the largest phi')

% plot the weighted average of complexes
numComplexes = length(phis_complex)-1; % Discard the trivial complex, i.e., the entire system.
type_of_heatmap = 1;
bipolar = 1;
[WeightedRatio_Complexes, AveragedPhi_Complexes] = AverageTopSubsets( complexes(1:end-1), phis_complex(1:end-1), N, numComplexes ); %Calculate average
figure(2)
make_ECoG_HeatMap( 'Chibi', target_ch, WeightedRatio_Complexes, type_of_heatmap, bipolar )
title(['Average of complexes. The number of complexes is ', num2str(numComplexes)])

% plot the weighted average of subsets
numTops = 10;
type_of_heatmap = 1;
bipolar = 1;
[WeightedRatio, AveragedPhi] = AverageTopSubsets( bipartitions2indices(Res.Z), Res.phi, N, numTops ); %Calculate average
figure(3)
make_ECoG_HeatMap( 'Chibi', target_ch, WeightedRatio, type_of_heatmap, bipolar )
title(['Average of subsets with the top', num2str(numTops)])

% Subsets and their phis
figure(4)
hoge = sortrows([Res.phi, Res.Z], -1);
subplot(2,1,1), imagesc(hoge(:,2:end)'),title('Subsets')
subplot(2,1,2), plot(hoge(:,1)), xlim([0.5 length(Res.phi)+0.5]),title('\Phi')