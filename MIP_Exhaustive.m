function [Z_MIP, phi_MIP] = MIP_Exhaustive( type_of_phi, Cov_X, Cov_XY, Cov_Y )
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

addpath(genpath('nextvector'))

%type_of_phi_str = string(type_of_phi);
N = size(Cov_X, 1);

IndexOutput_exhaustive = 1;
phi_MIP = Inf;
phis = zeros((2^N-2)/2,1);
idx = 1;
for iN = 1:floor(N/2)
    subcluster = 1:iN;
    if 2*iN == N
        nComb = nchoosek(N, iN)/2;
    else
        nComb = nchoosek(N, iN);
    end
    for j = 1:nComb
        Z = 2*ones(1, N);
        Z(subcluster) = 1;
        if strcmp( type_of_phi, 'MI1' )
            phi = phi_Gauss( type_of_phi, Z, Cov_X );
        else
            phi = phi_Gauss( type_of_phi, Z, Cov_X, Cov_XY, Cov_Y);
        end
        phis(idx) = phi;
        if phi < phi_MIP
            phi_MIP = phi;
            IndexOutput_exhaustive = subcluster;
        end
        subcluster = nextchoose(subcluster, N);
        idx = idx+1;
    end
end

Z_MIP = ones(1, N);
Z_MIP(IndexOutput_exhaustive) = 2;

end