%% Extract complexes from a toy network (Fig. 7 in Kitazono et al., 2022).
% If you set type_of_mat = 1, bidirectionality will be considered and 
% if you set type_of_mat = 2, bidirectionality will be ignored.

addpath(genpath('../../BidirectionallyConnectedCores'))
load('ToyNetworkData.mat')

type_of_mat = 1; % 1: Bidirectionality is considered, 2: Bidirectionality is ignored

if type_of_mat == 1
    W = ConnectionMatrix;
    G = digraph(W);
elseif type_of_mat == 2
    W = (ConnectionMatrix + ConnectionMatrix')/2;
    G = graph(W);
end


%% %%%%%%%%%%%%%%%%%%%%%%% Extract Complexes %%%%%%%%%%%%%%%%%%%%%%%%%
% Extract complexes
[complexes, w_mc_complexes, main_complexes, w_mc_main_complexes, Res] = HPC( W );

% Sort indices accotding to the hierarchical structure of complexes (See Fig. 6 in Kitazono et al., 2022)
[~, indices_sorted] = sortIndicesAccordingToHierarchicalStructure(Res);



%% %%%%%%%%%%%%%%%%%%%%%%% Visualize results %%%%%%%%%%%%%%%%%%%%%%%%%
%  Set figure parameters
mainSize = 0.80;
subSize = 0.025;
margin = (1-mainSize-subSize)/2;
gap_between_MainFigAndColorbar = 0.01;

Position_mat = [margin margin mainSize mainSize];
Position_clb = [margin+mainSize+gap_between_MainFigAndColorbar  margin subSize mainSize];

%% Visualize with network diagrams
% Network diagrams
figure('Position', [50 350 250 200]);
plot(G, 'XData', Xs, 'YData', Ys, 'NodeLabel', nodeNames, 'MarkerSize', 10, 'LineWidth', 2);

% Complexes
figure('Position', [300 350 250 200]);
colormap(parula)
plotComplexes(complexes, w_mc_complexes, '2D', 'Graph', G, 'XData', ...
Xs, 'YData', Ys, 'NodeLabel', nodeNames, 'MarkerSize', 10, 'LineWidth', 2);

figure('Position', [300 50 250 200]);
colormap(parula)
plotComplexes(complexes, w_mc_complexes, 'BirdsEye', 'Graph', G, 'XData', ...
Xs, 'YData', Ys, 'NodeLabel', nodeNames, 'MarkerSize', 10, 'LineWidth', 2);
zlabel('w^{mc}')


%% Visualize with connection matrices
% Connection matrix
figure('Position', [550 350 250 250]);
ax0 = axes('Position', Position_mat);
imagesc(ax0, W);
ax0.XTick = 1:10; ax0.YTick = 1:10;
xticklabels(nodeNames), yticklabels(nodeNames)
colormap(1-gray);
colorbar(ax0, 'Position', Position_clb);
title('Connection matrix')


% Visualize complexes by square areas.
figure('Position', [800 350 250 250]);
ax1 = axes('Position', Position_mat);
ComplexDiagram(ax1, complexes, w_mc_complexes, Res);
ax1.XTick = 1:10; ax1.YTick = 1:10;
xticklabels(nodeNames(indices_sorted)), yticklabels(nodeNames(indices_sorted))
colorbar(ax1, 'Position', Position_clb);
title('Sorted connection matrix')
% Superimpose the sorted connection matrix.
ax2 = axes('Position', Position_mat);
c = imagesc(ax2, W(indices_sorted, indices_sorted));
AlphaData = (W(indices_sorted, indices_sorted)~=0);
c.AlphaData = AlphaData;
colormap(ax2, 1-gray);
ax2.Visible = 'off';
% The order of the rows and columns may be partially different from Fig. 6,
% but this is not a problem, i.e., the extracted complexes are the same. 
% The order of the rows and columns is determined based on the hierarchical
% partitioning. In this hierarchical partitioning, if there are multiple
% min-cuts in a subnetwork, there is some arbitrariness in the order of
% partitionings. Accordingly, the order of rows and columns can be
% partially different. However, this arbitrariness does not affect the
% extracted complexes. 


%% Coreness
kv = coreness(complexes, w_mc_complexes, length(W));
figure('Position', [1050 350 200 150]);
plot(kv(indices_sorted), 'o')
title('Coreness')
xticks(1:length(W))
xticklabels(nodeNames(indices_sorted))
xlim([0.5 10.5])
ylim([-0.4 4.4])

