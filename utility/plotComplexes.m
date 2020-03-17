function plotComplexes(varargin)

[cax,args] = axescheck(varargin{:});
complexes = args{1};
phis = args{2};
plotType = args{3};
XData = args{4};
YData = args{5};
switch plotType
    case '3D'
        ZData = args{6};
        if length(args) > 6
            args = args(7:end);
        else
            args = {};
        end
    otherwise
        if length(args) > 5
            args = args(6:end);
        else
            args = {};
        end
end
nNodes = length(XData);

cLim = [];
for iArgs = 1:2:length(args)
    if strcmpi(args{iArgs}, 'CLim')
        cLim = args{iArgs+1};
        args = args([1:iArgs-1, iArgs+2:end]);
        break
    end
end
if isempty(args) || ~any(strcmpi('EdgeAlpha', args))
    args = [args, {'EdgeAlpha'}, {1}];
end
if isempty(args) || ~any(strcmpi('NodeLabel', args))
    args = [args, {'NodeLabel'}, {1:nNodes}];
end


nComplexes = length(phis);
[phis_sorted, index_phis_sorted] = sort(phis, 'ascend');
complexes_sorted = complexes(index_phis_sorted);

if isempty(cLim)
    cLim = [phis_sorted(1), phis_sorted(end)];
end

phis_sorted_rescaled = my_rescale(phis_sorted, cLim(1), cLim(2));
cmap = colormap;
nBinsColormap = size(cmap,1);
colorIndices = 1 + floor( (nBinsColormap-1)*phis_sorted_rescaled );
colors = cmap(colorIndices,:);
caxis(cLim)

hold on
G = ones(nNodes);
G(1:(nNodes+1):end) = 0;
G = graph(G);
for iComplexes = 1:nComplexes
    indices_temp = complexes_sorted{iComplexes};
    
    phi_temp = phis_sorted(iComplexes);
    G_temp = subgraph(G, indices_temp);
    subGPlotProps_temp = getSubGraphPlotProperties(G, indices_temp, args{:});
    
    switch plotType
        case '2D'
            ZData_temp = zeros(length(indices_temp), 1);
        case 'BirdsEye'
            ZData_temp = repelem(phi_temp, length(indices_temp));
        case '3D'
            ZData_temp = ZData(indices_temp);
    end
    plot(G_temp, 'NodeColor', colors(iComplexes,:), ...
        'EdgeColor', colors(iComplexes,:), ...
        'XData', XData(indices_temp), 'YData', YData(indices_temp), ...
        'ZData', ZData_temp, subGPlotProps_temp{:} );
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
for iIndices = 1:(length(indices)-1)
    idxOut = findedge(G, indices(iIndices), indices((iIndices+1):end));
    subGEdgeIndices = [subGEdgeIndices; nonzeros(idxOut)];
end

% [sOut,tOut] = findedge(G);
% nEdges = length(sOut); % nEdges = numedges(G);
% isMembSubG = zeros(nEdges, 1);
% for iEdges = 1:nEdges
%     if any(sOut(iEdges)==indices) && any(tOut(iEdges)==indices)
%         isMembSubG
%     end
% end

%subG = subgraph(G, indices);
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