function [EdgePhis, EdgeColors] = ComplexGraphs(X, Y, indices, phis, type_of_colormap, LineWidth, g, isshown, clim)
% Visualize complexes as superimposed graphs
%
% INPUT:
%    X, Y: X and Y coordinates of nodes
%    indices: The indices of complexes. (#complexes-by-1 cell array. Each cell contains indices of a subset.)
%    phis: The amount of integrated information for complexes. (#complexes-by-1 vector)
%    type_of_colormap: ex. 'parula'
%    LineWidth: max. width of edges
%    g: relative width of edges
%    isshown: show the image or not
%    clim: clim of colormap
%
% OUTPUT:
%    EdgePhis: the amount of integrated information of the complex where
%    each edge is included
%    EdgeColors: the color of each edge
%
% Jun Kitazono, 2018

nBinsColormap = 1024;
nElems = length(X);

nComplexes = length(phis);

[Phis_sorted, index_Phis_sorted] = sort(phis, 'ascend');
Complexes_sorted = indices(index_Phis_sorted);

if isempty(g) || max(g(:))<=0
    g = ones(nElems);
else
    max_g = max(g(:));
    g = g/max_g;
end


if nargin < 9
    Phis_sorted_rescaled = rescale(Phis_sorted);
    clim = [Phis_sorted(1), Phis_sorted(end)];
else
    if clim(1) ~= clim(2)
        Phis_sorted_rescaled = (Phis_sorted-clim(1))./(clim(2)-clim(1));
    else
        Phis_sorted_rescaled = 0;
    end
end

Colors = eval([type_of_colormap, '(nBinsColormap);']);
ColorIndices = 1 + floor( (nBinsColormap-1)*Phis_sorted_rescaled );
Colors = Colors(ColorIndices,:);

EdgePhis = zeros(nElems, nElems);
EdgeColors = cell(nElems, nElems);
for iComplexes = 1:nComplexes
    Pairs = nchoosek(Complexes_sorted{iComplexes}, 2);
    for iPairs = 1:size(Pairs, 1)
        Pair = Pairs(iPairs, :);
        EdgePhis(Pair(1), Pair(2)) = Phis_sorted(iComplexes);
        EdgeColors{Pair(1), Pair(2)} = Colors(iComplexes,:);
    end
end
EdgePhis = EdgePhis + EdgePhis';
%EdgeColors
%EdgeColors = EdgeColors + EdgeColors';

if isshown
    hold on
    Pairs = nchoosek(1:nElems, 2);
    EdgePhis_vec = zeros(size(Pairs, 1), 1);
    for iPairs = 1:size(Pairs, 1)
        Pair = Pairs(iPairs, :);
        EdgePhis_vec(iPairs) = EdgePhis(Pair(1), Pair(2));
    end
    [EdgePhis_vec_sorted, ind_sort] = sort(EdgePhis_vec, 'ascend');
    Pairs_sorted = Pairs(ind_sort,:);

    scatter(X, Y, 'Visible', 'off');
    for iPairs = 1:size(Pairs, 1)
        Pair = Pairs_sorted(iPairs,:);
        Xs = X(Pair);
        Ys = Y(Pair);
        if ~isempty(EdgeColors{Pair(1), Pair(2)}) && LineWidth*g(Pair(1), Pair(2))>0
            line(Xs, Ys, 'Color', EdgeColors{Pair(1), Pair(2)}, 'LineWidth', LineWidth*g(Pair(1), Pair(2)) );
            
            scatter(Xs, Ys, MarkerSize, EdgeColors{Pair(1), Pair(2)}, 'filled')
        end
    end
    
%     for i = 1:nElems
%         for j = i+1:nElems
%             
%             Xs = [X(i), X(j)];
%             Ys = [Y(i), Y(j)];
%             
%             line(Xs, Ys, 'Color', EdgeColors{i,j}, 'LineWidth', LineWidth);
%             
%         end
%     end
    
    % for iComplexes = 1:nComplexes
    %
    %     Pairs = nchoosek(Complexes_sorted{iComplexes}, 2);
    %     Color = Colors(iComplexes,:);
    % %     LineWidth = LineWidths(iComplexes);
    %     for iPairs = 1:size(Pairs,1)
    %         Pair = Pairs(iPairs, :);
    %         Xs = X(Pair);%[X(Pair(:,1)), X(Pair(:,2))];
    %         Ys = Y(Pair);%[Y(Pair(:,1)), Y(Pair(:,2))];
    %
    %         line(Xs, Ys, 'Color', Color, 'LineWidth', LineWidth);
    %     end
    %
    % end
    
    colormap(type_of_colormap);
    
    if clim(1)~=clim(2)
        caxis(clim)
        colorbar('Ticks', clim, 'TickLabels', {num2str(clim(1)), num2str(clim(end))})
    end
end

end

