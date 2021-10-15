function [complexes, w_mc_complexes, main_complexes, w_mc_main_complexes, Res] = HPC( W )
% Find complexes using Hierarchical Partitioning for Complex Search (HPC)
% Kitazono et al., 2021, bioRxiv. https://doi.org/10.1101/2021.07.12.452022
%
% INPUTS:
%           W: A connection matrix. The (i,j)-element indicate the weight
%              of the edge from the node i to the node j.
%
% OUTPUTS:
%    complexes: indices of elements in complexes (#complexes by 1 cell
%               array. Each cell contains indices of a complex).
%    w_mc_complexes: min-cut weights of complexes (#subsets by 1 vector).
%    main_complexes: indices of elements in main complexes
%    w_mc_main_complexes: min-cut weights of main complexes
%   
%    Res.Z: Each row of Z represents a complex candidate. This is a subset
%           belonging to \mathcal{V} in Kitazono et al., 2021. The i-th
%           element of a row has a nonzero value (1 or 2) when the element
%           is included in the candidate and zero when the element is NOT
%           included. The value 1 and 2 indicate the min-cut of a candidate. 
%           Ex.: Res.Z(i,:): [0, 0, 2, 1, 2, 0, 2]
%                             1  2  3  4  5  6  7
%                (The 3rd, 4th, 5th, and 7th elements are included in the
%                i-th candidate. Its min-cut divides it into {4} and {3,5,7}.) 
%           
%    Res.w_mc: min-cut weight of candidate complexes. The i-th element of 
%              Res.w_mc is min-cut weight of i-th candidate.
%
%    Res.parent: The parent of each candidate complex. When the i-th
%                element of Res.parent is j (i<j), this means that the j-th
%                candidate is the parent of i-th candidate in the
%                hierarchy. When i-th element is zero, this means the i-th
%                candidate is the root node (the whole system). 
%                
%    Res.isComplex: This indicates whether each candidate is a complex or
%                   not (1: a complex, 0: not a complex). 
%    Res.isMainComplex: This indicates whether each candidate is a main
%                       complex or not (1: a main complex, 0: not a main complex). 
%
% Jun Kitazono, 2021



%'HaoOrlin': an algorithm for directed graphs, which is used in 
%            Kitazono et al., 2021, bioRxiv.
%            https://doi.org/10.1101/2021.07.12.452022
%'NagamochiIbaraki': an algorithm for undirected graphs
if issymmetric( W )
    minCutAlgo = 'NagamochiIbaraki';
else
    minCutAlgo = 'HaoOrlin';
end


Res = Complex_RecursiveFunction( W, minCutAlgo );
[complexes, w_mc_complexes, isComplex, ...
    main_complexes, w_mc_main_complexes, isMainComplex] = find_Complexes_fromRes(Res);
Res.isComplex = isComplex;
Res.isMainComplex = isMainComplex;

end

function Res = Complex_RecursiveFunction( W, minCutAlgo )

N = length(W);

if N < 2
    Res.Z = [];
    Res.w_mc = [];
    Res.parent = 0;
else
    
    switch minCutAlgo
        case 'HaoOrlin'
            [Z_mc, w_mc] = hao_orlin_mincut( W );
        case 'NagamochiIbaraki'
            [Z_mc, w_mc] = nagamochi_ibaraki( W );
    end
    
    indices_L = find(Z_mc==1);
    indices_R = find(Z_mc==2);
    
    g_L = W(indices_L, indices_L);
    g_R = W(indices_R, indices_R);

    Res_L = Complex_RecursiveFunction( g_L, minCutAlgo );
    Res_R = Complex_RecursiveFunction( g_R, minCutAlgo );
    
    Res = Concatenate_Res(N, Z_mc, w_mc, indices_L, indices_R, Res_L, Res_R);
    
end

end

function Res = Concatenate_Res(N, Z, w_mc, indices_L, indices_R, Res_L, Res_R)

[nRows_L, ~] = size(Res_L.Z);
[nRows_R, ~] = size(Res_R.Z);
Res.Z = zeros(nRows_L+nRows_R+1, N);
Res.Z(1:nRows_R, indices_R) = Res_R.Z;
Res.Z((nRows_R+1):(nRows_R+nRows_L), indices_L) = Res_L.Z;
Res.Z(nRows_R+nRows_L+1,:) = Z;

Res.w_mc = [Res_R.w_mc; Res_L.w_mc; w_mc];

Res.parent = zeros(nRows_L+nRows_R+1, 1);
Res.parent(1:nRows_R, 1) = Res_R.parent;
Res.parent((nRows_R+1):(nRows_R+nRows_L),1) = nRows_R + Res_L.parent;

if nRows_R > 0
    Res.parent(nRows_R, 1) = nRows_L+nRows_R+1;
end
if nRows_R+nRows_L>0
    Res.parent(nRows_R+nRows_L, 1) = nRows_L+nRows_R+1;
end

end

function [complexes, w_mc_complexes, isComplex, main_complexes, w_mc_main_complexes, isMainComplex] = find_Complexes_fromRes(Res)

nSubsets = length(Res.w_mc);
w_mc_temp_max = zeros(nSubsets, 1);
isComplex = false(nSubsets, 1);
isMainComplex = false(nSubsets, 1);
mrac = zeros(nSubsets, 1);

w_mc_temp_max(nSubsets) = Res.w_mc(nSubsets);
isComplex(nSubsets) = true;
isMainComplex(nSubsets) = true;
mrac(nSubsets) = nSubsets;
for i = nSubsets-1: -1 :1
    if Res.w_mc(i) > w_mc_temp_max(Res.parent(i))
        w_mc_temp_max(i) = Res.w_mc(i);
        isComplex(i) = true;
        
        isMainComplex( mrac(Res.parent(i)) ) = false;
        isMainComplex(i) = true;
        mrac(i) = i;
    else
        w_mc_temp_max(i) = w_mc_temp_max(Res.parent(i));
        
        mrac(i) = mrac(Res.parent(i));
    end
end

[complexes, w_mc_complexes] = sort_Complexes(isComplex, Res);
[main_complexes, w_mc_main_complexes] = sort_Complexes(isMainComplex, Res);

end

function [complexes, w_mc] = sort_Complexes(isComplex, Res)

nComplexes = nnz(isComplex);
complexes = cell(nComplexes, 1);
Zs = Res.Z(isComplex,:);
for i = 1:nComplexes
    complexes{i} = find(Zs(i,:));
end

w_mc = Res.w_mc(isComplex);

[w_mc_sorted, idx_w_mc_sorted] = sort(w_mc, 'descend');
w_mc = w_mc_sorted;
complexes = complexes(idx_w_mc_sorted);

end