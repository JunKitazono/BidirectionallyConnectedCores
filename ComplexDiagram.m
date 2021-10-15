function [cax, indices_sorted] = ComplexDiagram(varargin)
% Ex. ComplexDiagram(complexes, w_mc_complexes, Res)


[cax,args] = axescheck(varargin{:});

complexes = args{1};
phis = args{2};
Res = args{3};

if isempty(cax) || ishghandle(cax,'axes')
    cax = newplot(cax);
    parent = cax;
else
    parent = cax;
    cax = ancestor(cax,'axes');
end


[phis_sorted, index_phis_sorted] = sort(phis, 'ascend');
complexes_sorted = complexes(index_phis_sorted);

%%
[~, indices_sorted] = sortIndicesAccordingToHierarchicalStructure(Res);
[~, order_inverse] = sort(indices_sorted);

nCs = length(phis);

nElems = size(Res.Z, 2);

%cla reset

complexBoxes = zeros(nElems);
for iCs = 1:nCs
    indices = complexes_sorted{iCs};
    indices_inverse = order_inverse(indices);
    minIndex = min(indices_inverse);
    maxIndex = max(indices_inverse);
    
    complexBoxes(minIndex:maxIndex, minIndex:maxIndex) = phis_sorted(iCs);
end
imagesc(cax, complexBoxes)
hold on


for iCs = 1:nCs
    indices = complexes_sorted{iCs};
    indices_inverse = order_inverse(indices);
    minIndex = min(indices_inverse);
    maxIndex = max(indices_inverse);

    minCoor = minIndex-0.5;
    maxCoor = maxIndex+0.5;
    
    xs_c = [minCoor maxCoor maxCoor minCoor minCoor];
    ys_c = [minCoor minCoor maxCoor maxCoor minCoor];
    line(cax, xs_c, ys_c, 'color', 'w')%colors(iCs,:))
end

end