% addpath(genpath('../PhiToolbox'))
% 
% nElems = 50;
% nElems_Sub = 10;
% Z = floor(((1:nElems)-1)/nElems_Sub)+1;
% w_intra = 1/nElems;
% w_inter = w_intra;
% ratios = 0:0.02:1;
% 
% As = make_ModularNetworkConnectivityMatrix(Z, w_intra, w_inter, ratios);
% 
% Cov_E = 0.01*eye(nElems);
% 
% Cov_Xs = cell(length(ratios), 1);
% Cov_XYs = cell(length(ratios), 1);
% Cov_Ys = cell(length(ratios), 1);
% probss = cell(length(ratios), 1);
% for i = 1:length(ratios)
%     [ probss{i}.Cov_X, probss{i}.Cov_XY, probss{i}.Cov_Y ] = symAtoCovs( As{i}, Cov_E );
%     probss{i}.number_of_elements = nElems;
% end
% 
% options.type_of_dist = 'Gauss';
% options.type_of_phi = 'MI1';
% options.type_of_MIPsearch = 'Queyranne';
% 
% toolboxes = ver;
% for i = 1:length(toolboxes)
%     if strcmp(toolboxes(i).Name, 'Parallel Computing Toolbox')
%         delete(gcp('nocreate'))
%         numCores = feature('numCores');
%         parpool(numCores);
%     end
% end
% 
% length_ratios = length(ratios);
% Complexes = cell(length(ratios), 1);
% phis = cell(length(ratios), 1);
% Res = cell(length(ratios), 1);
% main_Complexes = cell(length(ratios), 1);
% main_phis = cell(length(ratios), 1);
% 
% parfor i = 1:length_ratios
%     [Complexes{i}, phis{i}, Res{i}, main_Complexes{i}, main_phis{i}] = ...
%         Complex_Recursive_probs( probss{i}, options );
% end


inds = [2 20];
for i = 1:length(inds)
    Colors = hsv(length(main_Complexes{inds(i)}));
    NodeColors = 0.9*ones(nElems,3);
    for j = 1:length(main_Complexes{inds(i)})
        for k = 1:length(main_Complexes{inds(i)}{j})
            NodeColors(main_Complexes{inds(i)}{j}(k),:) = Colors(j,:);
        end
    end
    subplot(2,2,i), imagesc(As{inds(i)});
    subplot(2,2,2+i), plot(Gs{inds(i)}, 'NodeColor', NodeColors, 'MarkerSize', 10, 'NodeLabel', {});
    xticks([])
    yticks([])
end



