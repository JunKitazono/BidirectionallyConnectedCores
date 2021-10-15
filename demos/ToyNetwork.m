%% Extract complexes from a toy network (Fig. 7 in Kitazono et al., 2021).
% If you set type_of_mat = 1, bidirectionality will be considered and 
% if you set type_of_mat = 2, bidirectionality will be ignored.

addpath(genpath('../../BidirectionallyConnectedCores'))
load('ToyNetworkData.mat')

type_of_mat = 1; % 1: Bidirectionality is considered, 2: Bidirectionality is ignored

if type_of_mat == 1
    W = ConnectionMatrix;
elseif type_of_mat == 2
    W = (ConnectionMatrix + ConnectionMatrix')/2;
end


%% Extract complexes
% extract complexes
[complexes, w_mc_complexes, main_complexes, w_mc_main_complexes, Res] = HPC( W );

% Sort indices accotding to the hierarchical structure of complexes (See Fig. 6 in Kitazono et al., 2021)
[~, indices_sorted] = sortIndicesAccordingToHierarchicalStructure(Res);



%%  Set figure parameters
figPosition = 100+[0 0 300 300];
mainSize = 0.80;
subSize = 0.025;
margin = (1-mainSize-subSize)/2;
gap_between_MainFigAndColorbar = 0.01;

Position_mat = [margin margin mainSize mainSize];
Position_clb = [margin+mainSize+gap_between_MainFigAndColorbar  margin subSize mainSize];

%% Show the connection matrix
figure
set(gcf, 'Position', figPosition);
ax0 = axes('Position', Position_mat);
imagesc(ax0, W);
ax0.XTick = 1:10; ax0.YTick = 1:10;
xticklabels(nodeNames), yticklabels(nodeNames)
colormap(1-gray);
colorbar(ax0, 'Position', Position_clb);


%% Visualize results
% Visualize complexes by square areas.
figure
set(gcf, 'Position', figPosition);
ax1 = axes('Position', Position_mat);
ComplexDiagram(ax1, complexes, w_mc_complexes, Res);
ax1.XTick = 1:10; ax1.YTick = 1:10;
xticklabels(nodeNames(indices_sorted)), yticklabels(nodeNames(indices_sorted))
colorbar(ax1, 'Position', Position_clb);

% Superimpose the sorted connection matrix.
ax2 = axes('Position', Position_mat);
c = imagesc(ax2, W(indices_sorted, indices_sorted));
AlphaData = (W(indices_sorted, indices_sorted)~=0);
c.AlphaData = AlphaData;
colormap(ax2, 1-gray);
ax2.Visible = 'off';
% The order of the rows and columns may be partially different from Fig. 7,
% but this is not a problem (the extracted complexes are the same). 
% The order of the rows and columns is determined based on the hierarchical
% partitioning. In this hierarchical partitioning, if there are multiple
% min-cuts in a subnetwork, there is some arbitrariness in the order of
% partitionings. Accordingly, the order of rows and columns can be
% partially different. However, this arbitrariness does not affect the
% extracted complexes. 


% Coreness
kv = coreness(complexes, w_mc_complexes, length(W));
figure
plot(kv(indices_sorted), 'o')
title('Coreness')
xticks(1:length(W))
xticklabels(nodeNames(indices_sorted))
xlim([0.5 10.5])
ylim([-0.4 4.4])

