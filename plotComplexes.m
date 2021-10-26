function p = plotComplexes(varargin)
%% PLOTCOMPLEXES plots complexes as graphs.
% 
% PLOTCOMPLEXES(COMPLEXES,W_MC_COMPLEXES,PLOTTYPE,...) plots complexes as
% graphs. The color of each graph indicates the min-cut weight. The colors
% of graphs are determined by mapping the values in W_MC_COMPLEXES to the
% colors in the current colormap.
%   COMPLEXES is a cell array each cell of which contains the indices of
%      elements in a complex.
%   W_MC_COMPLEXES is a double vector with the same length as complexes
%      each cell of which indicates the min-cut weight of the complex. 
%   PLOTTYPE specifies a plot type. There are three options for PLOTTYPE, '2D', 'BirdsEye' and '3D'. 
%      '2D':PLOTCOMPLEXES(COMPLEXES,W_MC_COMPLEXES,'2D',...) plots
%           complexes as graphs in 2D space. The graphs are superimposed in
%           the ascending order of the min-cut weights.
%
%           Example:
%             PLOTCOMPLEXES(COMPLEXES,W_MC_COMPLEXES,'2D','XData',X,'YData',Y,...)
%             plots complexes as graphs by placing nodes at the coordinates
%             specified by the vectors X and Y. 
%
%             PLOTCOMPLEXES(COMPLEXES,W_MC_COMPLEXES,'2D','N',N,...) plots
%             complexes as graphs with a circular layout.  
%
%             PLOTCOMPLEXES(COMPLEXES,W_MC_COMPLEXES,'2D','Graph',G,...)
%             plots complexes as graphs with the same layout as that of
%             plot(G,...). 
%
%      'BirdsEye':PLOTCOMPLEXES(COMPLEXES,W_MC_COMPLEXES,'BirdsEye',...)
%                 plots complexes as graphs with a bird's eye view. Each
%                 graph is placed at the height specified by w_mc. The (X,Y)
%                 coordinates of the nodes are specified by the same way as
%                 in the '2D' case. 
%
%           Example:
%             PLOTCOMPLEXES(COMPLEXES,W_MC_COMPLEXES,'BidsEye','XData',X,'YData',Y,...)
%             PLOTCOMPLEXES(COMPLEXES,W_MC_COMPLEXES,'BidsEye','N',N,...)
%             PLOTCOMPLEXES(COMPLEXES,W_MC_COMPLEXES,'BidsEye','Graph',G,...)
%
%      '3D':PLOTCOMPLEXES(COMPLEXES,W_MC_COMPLEXES,'3D',...) plots
%           complexes as graphs in 3D space. The graphs are superimposed in
%           the ascending order of the min-cut weights.
%
%           Example:
%             PLOTCOMPLEXES(COMPLEXES,W_MC_COMPLEXES,'3D','XData',X,'YData',Y,'ZData',Z,...)
%             plots complexes as graphs by placing nodes at the coordinates
%             specified by the vectors X, Y and Z.
%
%
% PLOTCOMPLEXES(COMPLEXES,W_MC_COMPLEXES,'CLim',CLIM,...) specifies color
% limits. CLIM is a two-element vector of the form [cmin cmax], where cmax
% is greater than cmin. Values in W_MC_COMPLEXES that are less than or
% equal to cmin map to the first color in the colormap. Values greater than
% or equal to cmax map to the last color in the colormap. Values between
% cmin and cmax linearly map to the colormap. 
%
%
% PLOTCOMPLEXES(AX,...) plots into the axes with handle AX.
%
%
% H = plotComplexes(...) also returns a GraphPlot object for '2D' and '3D',
% and a GraphPlot object array of length(complexes) for 'BirdsEye'. Use the
% methods and properties of this object to inspect and adjust the plotted
% graph. 
%
% 
% PLOTCOMPLEXES accept the same Name-Value pair arguments as plot(G), 
% except Linespec (e.g., '-.or'). 
% The previous syntaxes can also contain one or more Name-Value pair
% arguments that specify additional properties. 
%
% See also the link below.
% https://www.mathworks.com/help/matlab/ref/graph.plot.html
% https://www.mathworks.com/help/matlab/ref/matlab.graphics.chart.primitive.graphplot-properties.html
%
%
% Examples:
%     plotComplexes(complexes,w_mc_complexes,'2D', 'XData',X, 'YData',Y, 'LineWidth', 2, 'CLim', [0.1, 0.4])
%     plotComplexes(main_complexes,w_mc_main_complexes,'2D', 'XData',X, 'YData',Y)
%     plotComplexes(complexes, w_mc_complexes, 'BirdsEye', 'Graph', G, 'LineWidth', 2, 'Layout', 'force')
%     plotComplexes(complexes, w_mc_complexes, '3D','XData',X, 'YData',Y, 'ZData',Z, 'Marker','d','LineStyle','--')

% Jun Kitazono, 2020

[cax,args] = axescheck(varargin{:});
complexes = args{1};
w_mc = args{2};
plotType = args{3};
args = args(4:end);

if isempty(cax) || ishghandle(cax,'axes')
    cax = newplot(cax);
    parent = cax;
else
    parent = cax;
    cax = ancestor(cax,'axes');
end

[args, cLim, isCLim] = ExtractFromArgs(args, 'CLim');
[args, G, isGraph] = ExtractFromArgs(args, 'Graph');
[args, nNodes, isN] = ExtractFromArgs(args, 'N');
[args, XData, isX] = ExtractFromArgs(args, 'XData');
[args, YData, isY] = ExtractFromArgs(args, 'YData');
[args, ZData, isZ] = ExtractFromArgs(args, 'ZData');

% Check X and Y
if xor(isX, isY)
    error('You must use ''XData'' with ''YData''.')
end

% Check nNodes
if ~any([isN, isX, isGraph])
    error('Please specify ''N'', ''XData'', or ''Graph'' to determine the number of nodes.')
end
if isN && (isX || isGraph)
    error('You cannot use ''N'' with ''XData'' or ''Graph''.')
end
if ~isN && isX && isGraph
    if ~isequal(length(XData), numnodes(G))
        error('The length(XData) must be the same as numnodes(G).')
    end
end

% Check Layout
if any(strcmpi('Layout', args)) && ~isGraph
    error('You must use ''Layout'' with ''Graph''.')
end
if any(strcmpi('Layout', args)) && isX
    error('You cannot use ''Layout'' with ''X''.')
end

if ~isN
    if isGraph
        nNodes = numnodes(G);
    else % isX
        nNodes = length(XData);
    end
end

if isGraph
    G = rmedge(G, 1:nNodes, 1:nNodes); % Remove self loop for visualization
    adjG = adjacency(G);
else
    adjG = ones(nNodes);
    adjG(1:(nNodes+1):end) = 0;
end

switch plotType
    case {'2D', 'BirdsEye'}
        if ~(isX && isY)
            if isGraph
                f = figure;
                h = plot(G, args{:});
                XData = h.XData; YData = h.YData;
                close(f)
                clear h f
                args = ExtractFromArgs(args, 'Layout');
            elseif isN
                thetas = 2*pi*(1:nNodes)/nNodes;
                XData = cos(thetas); YData = sin(thetas);
            else
                error('Please specify (XData, YData), Graph, or N.') 
            end
        end
    case '3D'
        if ~(isX && isY && isZ)
            error('Please specify ''XData'', ''YData'', and ''ZData''.')
        end
end

args = InsertDefault2Args(args, 'EdgeAlpha', 1);
args = InsertDefault2Args(args, 'NodeLabel', (1:nNodes)');


nComplexes = length(w_mc);
[phis_sorted, index_phis_sorted] = sort(w_mc, 'ascend');
complexes_sorted = complexes(index_phis_sorted);

adjComps_all = zeros(nNodes);
indices_all = [];
for iComplexes = 1:nComplexes
    indices_temp = complexes_sorted{iComplexes};
    adjComps_all(indices_temp, indices_temp) = adjG(indices_temp, indices_temp);
    indices_all = union(indices_all, indices_temp);
end
if issymmetric(adjComps_all)
    GComps = graph(adjComps_all);
else
    GComps = digraph(adjComps_all);
end

edges_upper_union = [];
edges_upper = cell(nComplexes, 1);
for iComplexes = nComplexes:-1:1
    edges_upper{iComplexes} = edges_upper_union;
    
    indices_temp = complexes_sorted{iComplexes};
    [x,y] = meshgrid(indices_temp, indices_temp);
    x = x(:); y = y(:);
    edges_temp = nonzeros(findedge(GComps, x, y));
    edges_upper_union = union(edges_upper_union, edges_temp);
    
end

%% Setting colors
if ~isCLim
    cLim = [phis_sorted(1), phis_sorted(end)];
end

phis_sorted_rescaled = my_rescale(phis_sorted, cLim(1), cLim(2));
cmap = colormap;
nBinsColormap = size(cmap,1);
colorIndices = 1 + floor( (nBinsColormap-1)*phis_sorted_rescaled );
colors = cmap(colorIndices,:);
if cLim(1)~=cLim(2)
    caxis(cLim)
else
    caxis([cLim(1), cLim(1)+1])
end

%% Plot
hold on
switch plotType
    case {'2D'}
        [G_WithAttribute, GPlotProps_attributed, GPlotProps_remaining] = addGraphPlotProperties2G(GComps, args{:});
        
        for iComplexes = 1:nComplexes
            idx_temp = complexes_sorted{iComplexes};
            G_temp = rmedge(G_WithAttribute, edges_upper{iComplexes});
            G_temp = subgraph(G_temp, idx_temp);
            
            subGPlotProps_temp = getGraphPlotPropertiesFromG(G_temp, GPlotProps_attributed{:});
            
            hObj(iComplexes) = ...
                plot(cax, G_temp, 'Parent', parent,...
                'NodeColor', colors(iComplexes,:), ...
                'EdgeColor', colors(iComplexes,:), ...
                'XData',XData(idx_temp),...
                'YData',YData(idx_temp), ...
                subGPlotProps_temp{:}, GPlotProps_remaining{:} );
        end
    case {'BirdsEye'}
        for iComplexes = 1:nComplexes
            indices_temp = complexes_sorted{iComplexes};
            
            phi_temp = phis_sorted(iComplexes);
            G_temp = subgraph(GComps, indices_temp);
            subGPlotProps_temp = getSubGraphPlotProperties(GComps, indices_temp, args{:});
            ZData_temp = repelem(phi_temp, length(indices_temp));
            
            hObj(iComplexes) = ...
                plot(cax, G_temp, 'Parent', parent,...
                'NodeColor', colors(iComplexes,:), ...
                'EdgeColor', colors(iComplexes,:), ...
                'XData',XData(indices_temp),...
                'YData',YData(indices_temp), ...
                'ZData',ZData_temp, subGPlotProps_temp{:} );
        end
    case '3D'
        GAll = subgraph(GComps, indices_all);
        subGPlotProps_all = getSubGraphPlotProperties(GComps, indices_all, args{:});
        
        hObj = plot(cax, GAll, 'Parent', parent,...
            'NodeColor', 'k', 'EdgeColor', 'k', ...
            'XData',XData(indices_all), 'YData',YData(indices_all),...
            'ZData',ZData(indices_all), subGPlotProps_all{:});
        for iComplexes = 1:nComplexes
            indices_temp = complexes_sorted{iComplexes};
            G_temp = zeros(nNodes);
            G_temp(indices_temp, indices_temp) = adjG(indices_temp, indices_temp);
            if issymmetric(G_temp)
                G_temp = graph(G_temp);
            else
                G_temp = digraph(G_temp);
            end
            G_temp = subgraph(G_temp, indices_all);
            highlight(cax.Children(1), G_temp, 'NodeColor', colors(iComplexes,:), ...
                'EdgeColor', colors(iComplexes,:));
        end
end

switch plotType
    case {'3D', 'BirdsEye'}
        view(-10, 20)
end

hold off

if nargout > 0
    p = hObj;
end
        
end

function [args, val, isNamedArg] = ExtractFromArgs(args, name)

indVal = find(strcmpi(name, args), 1);
if ~isempty(indVal)
    val = args{indVal+1};
    args = args([1:indVal-1, indVal+2:end]);
    isNamedArg = true;
else
    val =[];
    isNamedArg = false;
end

end

function args = InsertDefault2Args(args, name, defaultVal)

if isempty(args) || ~any(strcmpi(name, args))
    args = [args, {name}, {defaultVal}];
end

end

function phis_rescaled = my_rescale(phis, minPhi, maxPhi)
phis_rescaled = zeros(size(phis));
if minPhi ~= maxPhi
    phis_rescaled = (phis-minPhi)./(maxPhi-minPhi);
end
phis_rescaled = min(phis_rescaled, 1);
phis_rescaled = max(phis_rescaled, 0);

end


function subGPlotProps = getSubGraphPlotProperties(varargin)

G = varargin{1};
indices = varargin{2};
GPlotProps = varargin(3:end);

subGEdgeIndices = [];
for iIndices = 1:length(indices)
    %idxOut = findedge(G, indices(iIndices), indices(iIndices:end));
    idxOut = findedge(G, indices(iIndices), indices);
    subGEdgeIndices = [subGEdgeIndices; nonzeros(idxOut)];
end
subGEdgeIndices = unique(subGEdgeIndices);

subGPlotProps = cell(size(GPlotProps));
for i = 1:2:length(GPlotProps)
    
    subGPlotProps{i} = GPlotProps{i};
    
    if size(GPlotProps{i+1}) > 1 % Color, ()*3 mat
        if startsWith(GPlotProps{i}, 'Node', 'IgnoreCase', true) % Nodes
            subGPlotProps{i+1} = GPlotProps{i+1}(indices,:);
        elseif startsWith(GPlotProps{i}, 'Edge', 'IgnoreCase', true) % Edges
            subGPlotProps{i+1} = GPlotProps{i+1}(subGEdgeIndices,:);
        end
    elseif length(GPlotProps{i+1}) > 1 && ~isa(GPlotProps{i+1}, 'char') % vec other than char or callbacks
        if startsWith(GPlotProps{i}, 'Node') || startsWith(GPlotProps{i}, 'Marker') % Nodes
            subGPlotProps{i+1} = GPlotProps{i+1}(indices);
        elseif startsWith(GPlotProps{i}, 'Edge', 'IgnoreCase', true) || ...
                startsWith(GPlotProps{i}, 'Line', 'IgnoreCase', true) || ...
                startsWith(GPlotProps{i}, 'Arrows', 'IgnoreCase', true) % Edges
            subGPlotProps{i+1} = GPlotProps{i+1}(subGEdgeIndices);
        end
    else
        subGPlotProps{i+1} = GPlotProps{i+1};
    end
        
end

end

function [G_WithAttribute, GPlotProps_attributed, GPlotProps_remaining] = addGraphPlotProperties2G(varargin)

G_WithAttribute = varargin{1};
GPlotProps = varargin(2:end);

GPlotProps_attributed = cell(0);
GPlotProps_remaining = cell(0);

for i = 1:2:length(GPlotProps)
    
    if size(GPlotProps{i+1}) > 1 % Color, ()*3 mat
        if startsWith(GPlotProps{i}, 'Node', 'IgnoreCase', true) % Nodes
            eval(['G_WithAttribute.Nodes.',GPlotProps{i}, ' = GPlotProps{i+1};'])
        elseif startsWith(GPlotProps{i}, 'Edge', 'IgnoreCase', true) % Edges
            eval(['G_WithAttribute.Edges.',GPlotProps{i}, ' = GPlotProps{i+1};'])
        end
        GPlotProps_attributed = [GPlotProps_attributed, GPlotProps(i), GPlotProps(i+1)];
    elseif length(GPlotProps{i+1}) > 1 && ~isa(GPlotProps{i+1}, 'char') % vec other than char or callbacks
        if startsWith(GPlotProps{i}, 'Node') || startsWith(GPlotProps{i}, 'Marker') % Nodes
            eval(['G_WithAttribute.Nodes.',GPlotProps{i}, ' = GPlotProps{i+1};'])
        elseif startsWith(GPlotProps{i}, 'Edge', 'IgnoreCase', true) || ...
                startsWith(GPlotProps{i}, 'Line', 'IgnoreCase', true) || ...
                startsWith(GPlotProps{i}, 'Arrows', 'IgnoreCase', true) % Edges
            eval(['G_WithAttribute.Edges.',GPlotProps{i}, ' = GPlotProps{i+1};'])
        end
        GPlotProps_attributed = [GPlotProps_attributed, GPlotProps(i), GPlotProps(i+1)];
    else
        GPlotProps_remaining = [GPlotProps_remaining, GPlotProps(i), GPlotProps(i+1)];
    end
end

end

function GPlotProps = getGraphPlotPropertiesFromG(varargin)

G_WithAttribute = varargin{1};
GPlotProps_attributed = varargin(2:end);

GPlotProps = GPlotProps_attributed;
for i = 1:2:length(GPlotProps_attributed)
    if size(GPlotProps{i+1}) > 1 % Color, ()*3 mat
        if startsWith(GPlotProps{i}, 'Node', 'IgnoreCase', true) % Nodes
            eval(['GPlotProps{i+1} = G_WithAttribute.Nodes.', GPlotProps{i}, ';'])
        elseif startsWith(GPlotProps{i}, 'Edge', 'IgnoreCase', true) % Edges
            eval(['GPlotProps{i+1} = G_WithAttribute.Edges.', GPlotProps{i}, ';'])
        end
    elseif length(GPlotProps{i+1}) > 1 && ~isa(GPlotProps{i+1}, 'char') % vec other than char or callbacks
        if startsWith(GPlotProps{i}, 'Node') || startsWith(GPlotProps{i}, 'Marker') % Nodes
            eval(['GPlotProps{i+1} = G_WithAttribute.Nodes.', GPlotProps{i}, ';'])
        elseif startsWith(GPlotProps{i}, 'Edge', 'IgnoreCase', true) || ...
                startsWith(GPlotProps{i}, 'Line', 'IgnoreCase', true) || ...
                startsWith(GPlotProps{i}, 'Arrows', 'IgnoreCase', true) % Edges
            eval(['GPlotProps{i+1} = G_WithAttribute.Edges.', GPlotProps{i}, ';'])
        end
    end
end

end