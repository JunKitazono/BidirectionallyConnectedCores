function ComplexDiagram(varargin)
%complexes, phis, Res, w

% if nargin < 1,error(message('nnet:Args:NotEnough'));end
% if nargin < 2, max_m = max(max(abs(w))); end
% if nargin < 3, min_m = max_m / 100; end
% if max_m == min_m, max_m = 1; min_m = 0; end
[cax,args] = axescheck(varargin{:});

complexes = args{1};
phis = args{2};
Res = args{3};
w = args{4};

if isempty(cax) || ishghandle(cax,'axes')
    cax = newplot(cax);
    parent = cax;
else
    parent = cax;
    cax = ancestor(cax,'axes');
end

% [args, boxLim, isBoxLim] = ExtractFromArgs(args, 'BoxLim');
% if ~isBoxLim
%     max_m = max(max(w));
%     min_m = min(min(w));
%     if max_m == min_m, max_m = 1; min_m = 0; end
%     boxLim = [min_m, max_m];
% end


[phis_sorted, index_phis_sorted] = sort(phis, 'ascend');
complexes_sorted = complexes(index_phis_sorted);

[args, pixelAlphaLim, isPixelAlphaLim] = ExtractFromArgs(args, 'pixelAlphaLim');
if ~isPixelAlphaLim
    pixelAlphaLim = [0, max(max(w))];
end
AlphaData = my_rescale(w, pixelAlphaLim(1), pixelAlphaLim(2));


%% Setting colors
[args, cLim, isCLim] = ExtractFromArgs(args, 'CLim');
if ~isCLim
    cLim = [phis_sorted(1), phis_sorted(end)];
end
if cLim(1)~=cLim(2)
    caxis(cLim)
else
    caxis([cLim(1), cLim(1)+1])
end


%%
[Z_sorted, indices_sorted] = sortIndicesWithHierarchicalStructure(Res, 0);
[~, order_inverse] = sort(indices_sorted);

nCs = length(phis);

[S,R] = size(w);

%cla reset

complexBoxes = zeros(S, R);
for iCs = 1:nCs
    indices = complexes_sorted{iCs};
    indices_inverse = order_inverse(indices);
    minIndex = min(indices_inverse);
    maxIndex = max(indices_inverse);
    
    complexBoxes(minIndex:maxIndex, minIndex:maxIndex) = phis_sorted(iCs);
end
imagesc(cax, complexBoxes)
hold on

% 
% % DEFINE POSITIVE BOX
% xn = [-1 -1 +1 +1 -1]*0.5;
% yn = [-1 +1 +1 -1 -1]*0.5;


%hold on
%%%%
% set(cax,'xlim',[0 R]+0.5);
% set(cax,'ylim',[0 S]+0.5);
% set(cax,'xlimmode','manual');
% set(cax,'ylimmode','manual');
% xticks = get(cax,'xtick');
% set(cax,'xtick',xticks(find(xticks == floor(xticks))))
% yticks = get(cax,'ytick');
% set(cax,'ytick',yticks(find(yticks == floor(yticks))))
% set(cax,'ydir','reverse');
%%%%
% if get(0,'screendepth') > 1
%   set(cax,'color',[1 1 1]*.5);
%   set(gcf,'color',[1 1 1]*.3);
% end

% for i=1:S
%     for j=1:R
%         % m = sqrt( (w_sorted(i,j)-min_m)/(max_m-min_m) );
%         % m = min(m,max_m)*0.95;
%         m = sqrt( my_rescale(w_sorted(i,j), min_m, max_m) ) * 0.95;
%         if m > 0
%             fill(xn*m+j,yn*m+i,[0 0 0])
%         end
%     end
% end


%colormap([1 1 1; 0 0 0])

isEdge = (w==0);
edgeC = uint8(repmat(isEdge*255, 1, 1, 3));
image(cax, edgeC, 'AlphaData', AlphaData, 'AlphaDataMapping', 'none')

%phis_sorted_rescaled = my_rescale(phis_sorted, cLim(1), cLim(2));
% cmap = colormap;
% nBinsColormap = size(cmap,1);
% colorIndices = 1 + floor( (nBinsColormap-1)*phis_sorted_rescaled );
% colors = cmap(colorIndices,:);
for iCs = 1:nCs
    indices = complexes_sorted{iCs};
    indices_inverse = order_inverse(indices);
    minIndex = min(indices_inverse);
    maxIndex = max(indices_inverse);
    %boxSize = maxIndex - minIndex + 1;
    
    minCoor = minIndex-0.5;
    maxCoor = maxIndex+0.5;
    
    xs_c = [minCoor maxCoor maxCoor minCoor minCoor];
    ys_c = [minCoor minCoor maxCoor maxCoor minCoor];
    line(cax, xs_c, ys_c, 'color', 'w')%colors(iCs,:))
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

function phis_rescaled = my_rescale(phis, minPhi, maxPhi)
phis_rescaled = zeros(size(phis));
if minPhi ~= maxPhi
    phis_rescaled = (phis-minPhi)./(maxPhi-minPhi);
end
phis_rescaled = min(phis_rescaled, 1);
phis_rescaled = max(phis_rescaled, 0);

end