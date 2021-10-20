function [Z_MIP, phi_MIP] = nagamochi_ibaraki(g)
% FUNCTION: nagamochi_ibaraki.m
% PURPOSE: Find the minimum cut in the undirected graph
% 
% INPUT: 
%       g       adjacent matrix of graph. gi,j is weight 
%               of the edge connecting nodes i,j.
%               the matrix must be symmetric.
% 
% OUTPUTS:
%       Z_MIP   the estimated MIP(grouped into 1 and 2)
%       
%       phi_MIP the value of mincut of the graph which corresponds to 
%               integrated information at the MIP
% 
% EXAMPLE:
%       test_data = [0 0.2 0.5; 0.2 0 0.1; 0.5 0.1 0]; %give an adjacent
%       matrix
%       [Z_MIP, phi_MIP] = nagamochi_ibaraki(test_data);
%       
%   result:
%       Z_MIP = [2 1 2], phi_MIP = 0.3000
%       means the graph was split to two groups: {1,3}, {2}
%       minimum cut weight value is 0.3000
%
% Complexity:
%       O(|V||E|log(|V|))
%
% USAGE:
%       this function uses open source library called LEMON
%       compile nagamochi_ibaraki_c.cpp file by following command at the
%       proper directory
%       mex BidirectionallyConnectedCores/tools/lemon/nagamochi_ibaraki_c.cpp -I/usr/local/include/ -L/usr/local/lib/ -lemon
%       You must install lemon library to your local environment
% 
% Yuma Aoki, 2020

if g ~= transpose(g)
    error('input matrix must be symmetric')
end
    

[Z_MIP, phi_MIP] = hao_orlin_mincut_c(g);
Z_MIP = Z_MIP + 1;


