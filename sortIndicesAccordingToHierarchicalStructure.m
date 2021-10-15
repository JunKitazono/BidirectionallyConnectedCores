function [Z_sorted, indices_sorted] = sortIndicesAccordingToHierarchicalStructure(Res)

if isfield(Res, 'phi')
    phi = Res.phi;
elseif isfield(Res, 'w_mc')
    phi = Res.w_mc;
end

Z = Res.Z;


indices_sorted = 1:size(Z, 2);
for i = 1:length(phi)
    Z_temp = Z(i,:);
    
    Z_nonzeros = nonzeros(Z_temp);%Z(i, Z(i,:)>0);
    [~, arrangement_temp] = sort(Z_nonzeros); % [B,I]=sort(A); B=A(I)
    
    indices_sorted_sub = indices_sorted(Z_temp>0);
    indices_sorted_sub = indices_sorted_sub(arrangement_temp);
    indices_sorted(Z_temp>0) = indices_sorted_sub;
    
end

Z_sorted = Z(:, indices_sorted);

end