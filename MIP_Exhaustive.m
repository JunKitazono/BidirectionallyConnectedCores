function [Z_MIP, phi_MIP] = MIP_Exhaustive( type_of_dist, type_of_phi, X, params, probs)
%-----------------------------------------------------------------------
% FUNCTION: MIP_Exhaustive.m
% PURPOSE: Find the minimum information partition (MIP) by Exhaustive Search
%
% INPUTS:   
%           type_of_phi:
%              'SI': phi_H, stochastic interaction
%              'Geo': phi_G, information geometry version
%              'star': phi_star, based on mismatched decoding
%              'MI': Multi (Mutual) information, I(X_1, Y_1; X_2, Y_2)
%              'MI1': Multi (Mutual) information. I(X_1; X_2). (IIT1.0)
%           Cov_X: covariance of data X (past, t-tau)
%           Cov_XY: cross-covariance of X (past, t-tau) and Y (present, t)
%           Cov_Y: covariance of data Y (present, t)
%              
%
% OUTPUT:
%           Z_MIP: the MIP
%           phi_MIP: the amount of integrated information at the MIP
%-----------------------------------------------------------------------

if nargin < 5
    probs = data_to_probs(type_of_dist,X,params);
end

N = size(X,1);
all_comb = power_set(2:N);
nComb = length(all_comb);
phis = zeros(nComb,1);

parfor i=1: nComb
    subcluster = all_comb{i};
    Z = ones(1,N);
    Z(subcluster) = 2;
    % compute phi
    phis(i) = phi_comp(type_of_dist, type_of_phi, Z, X, params, probs);
end

[phi_MIP, MIP_ind] = min(phis);

Z_MIP = ones(1,N);
Z_MIP(all_comb{MIP_ind}) = 2;

end
