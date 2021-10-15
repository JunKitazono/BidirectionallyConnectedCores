function kv = coreness( complexes, w_mc_complexes, nelems )
% Take maximum phi of each element. 
%
% INPUT:
%    complexes: indices of elements in complexes (#complexes by 1 cell
%    array. Each cell contains indices of a complex).
%    w_mc_complexes: min-cut weights of complexes (#complexes by 1 vector).
%    nelems: The number of elements in the entire network.
%
% OUTPUT:
%    kv: corenss (1 by nelems vector). The coreness of node v is kv if an
%    element is included in a complex with w_mc = kv but not included in
%    any complex with w_mc > kv.
%
% Jun Kitazono, 2021

[w_mc_sort, index_phis_sort] = sort(w_mc_complexes, 'descend');

kv = zeros(1, nelems);
for i = length(w_mc_complexes): -1: 1
    kv(1, complexes{index_phis_sort(i)}) = w_mc_sort(i);
end

end