%% Extract complexes from a mouse connectome (Fig. 8 in Kitazono et al., 2022)
% If you set type_of_mat = 1, bidirectionality will be considered and 
% if you set type_of_mat = 2, bidirectionality will be ignored.

addpath(genpath('../../BidirectionallyConnectedCores'))
load('MouseConnectome.mat')

type_of_mat = 1; % 1: Bidirectionality is considered, 2: Bidirectionality is ignored

if type_of_mat == 1
    W = ConnectionMatrix;
elseif type_of_mat == 2
    W = (ConnectionMatrix + ConnectionMatrix')/2;
end

%% Extract Complexes
% extract complexes
[complexes, w_mc_complexes, main_complexes, w_mc_main_complexes, Res] = HPC( W );

% Sort indices accotding to the hierarchical structure of complexes (See Fig. 6 in Kitazono et al., 2022)
[~, indices_sorted] = sortIndicesAccordingToHierarchicalStructure(Res);

kv = coreness(complexes, w_mc_complexes, size(W,1));

%% Draw figures

% Set figure parameters
figPosition = 100+[0 0 400 400];

mainSize = 0.80;
subSize = 0.025;
margin = (1-mainSize-subSize)/2;
mainMarginL = margin + subSize;
gap_between_MainFigAndColorbar = 0.01;

Position_mat = [mainMarginL  margin mainSize mainSize];
Position_clb = [mainMarginL+mainSize+gap_between_MainFigAndColorbar  margin subSize mainSize];
Position_clb_top = [mainMarginL margin+mainSize mainSize subSize];
Position_clb_left = [margin margin subSize mainSize];

clims = [0 10^(-0.5)]; % clim for edge weights

% 
for Original_or_Sorted = ["Original", "Sorted"]

switch Original_or_Sorted
    case "Original"
        indexOrder = 1:size(W,1);
    case "Sorted"
        indexOrder = indices_sorted;
end
 
figure
set(gcf, 'Position', figPosition)

    
% Draw the connection matrix
ax_W = axes('Position', Position_mat);
c = imagesc(ax_W, W(indexOrder, indexOrder), clims);
AlphaData = rescale(W(indexOrder, indexOrder), 'InputMin', clims(1), 'InputMax', clims(2));
c.AlphaData = AlphaData; % set transparency
colormap(ax_W, zeros(256,3))
ax_W.XTick  = [];
ax_W.YTick = [];


switch Original_or_Sorted
    case "Original"
        % Color bar for edge weights
        clb = colorbar(ax_W, 'Position', Position_clb);
        drawnow
        cdata = clb.Face.Texture.CData; % Set transparency of the colorbar for edge weights
        cdata(end,:) = uint8(0:255);
        clb.Face.Texture.ColorType = 'truecoloralpha';
        pause(1) % <- This line is necessary to enable the transparency setting for some unknown reason...
        clb.Face.Texture.CData = cdata;
    case "Sorted"
        ax_W.Visible = 'off';
                
        % Visualize the structure of complexes
        ax_c = axes('Position', Position_mat);
        ComplexDiagram(ax_c, complexes, w_mc_complexes, Res);
        ax_c.XTick  = [];
        ax_c.YTick = [];
        uistack(ax_c, 'bottom')
        
        % Color bar for w_mc
        clb = colorbar(ax_c, 'Position', Position_clb);
end


% Color bar that indicates major brain regions (the top one)
ax_clb_top = axes('Position', Position_clb_top);
image(ax_clb_top, fliplr(color_MajorRegionColorBar(indexOrder,:,:)))
ax_clb_top.XTick  = [];
ax_clb_top.YTick = [];
camroll(90)

% Color bar that indicates major brain regions (the left one) 
ax_clb_left = axes('Position', Position_clb_left); 
image(ax_clb_left, color_MajorRegionColorBar(indexOrder,:,:))
ax_clb_left.XTick  = [];
ax_clb_left.YTick = [];


end


% Violin plot of coreness values.
figure
%h_Violin = violinplot(kv, [convertCharsToStrings(cellstr(RegionInfo.MajorRegion)); repmat({'All'},213,1)], 'GroupOrder',[cellstr(MajorRegionList)', {'All'}]);

h_Violin = violinplot(kv, [cellstr(MajorRegionList(MajorRegionNumber(1:213))); repmat({'All'},213,1)], 'GroupOrder',[cellstr(MajorRegionList)', {'All'}]);
for i = 1:14
    if i <= 13
        h_Violin(i).ViolinColor = double(color_MajorRegionList(i, :, :))/255;
        h_Violin(i).EdgeColor = double(color_MajorRegionList(i, :, :))/255;
    else
        h_Violin(i).ViolinColor = [0.5 0.5 0.5];
        h_Violin(i).EdgeColor = [0.5 0.5 0.5];
    end
   h_Violin(i).ViolinAlpha = 0.2;
   h_Violin(i).ViolinPlot.EdgeAlpha = 0.2;   
   h_Violin(i).ScatterPlot.MarkerFaceAlpha = 0.8;
   
end
xtickangle(45)
