function [Z_MIP, phi_MIP] = hao_orlin_mincut(g)
% FUNCTION: hao_orlin_mincut.m
% PURPOSE: Find the minimum cut in the directed graph
% 
% INPUT: 
%       g       adjacent matrix of graph. gi,j is weight 
%               of the edge connecting nodes i,j.
% 
% OUTPUTS:
%       Z_MIP   the estimated MIP
%       
%       phi_MIP the value of mincut of the graph which corresponds to 
%               integrated information at the MIP
% 
% EXAMPLE:
%       test_data = [0 0.2 0.5; 0.2 0 0; 0.5 0.1 0]; %give an adjacent
%       matrix
%       [Z_MIP, phi_MIP] = hao_orlin_mincut(g);
%       
%   result:
%       Z_MIP = [2 1 2], phi_MIP = 0.2000
%       means the graph was split to two groups: {1,3}, {2}
%       minimum cut weight value is 0.2000
% Complexity:
%       O(|V|^2 * sqrt(|E|))
% USAGE:
%       this function uses open source library called LEMON
%       compile hao_orlin_mincut.cpp file by following command at the
%       proper directory
%       mex PhiToolbox/tools/lemon/hao_orlin_mincut_c.cpp -I/usr/local/include/ -L/usr/local/lib/ -lemon
%       You must install lemon library to your local environment
% 
% Yuma Aoki, 2019


[Z_MIP, phi_MIP] = hao_orlin_mincut_c(g);
Z_MIP = Z_MIP + 1;


