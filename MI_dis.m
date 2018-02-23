function MI = MI_dis(Z, N_st, p)

%-------------------------------------------------------------------------------------------------
% PURPOSE: calculate mutual (multi) information between subsystems in discretized data
%
% INPUTS:
%   Z: partition
%   N_st: number of states in each unit
%   p: probability distribution of time series X
%
% OUTPUTS:
%   MI: mutual information between subsystems 
%-------------------------------------------------------------------------------------------------
%
% Masafumi Oizumi, 2018

N = length(Z);
N_c = max(Z);
H_part = zeros(N_c,1);
for i=1: N_c
    sub_index = find(Z==i);
    q = marginalize(p,sub_index,N,N_st);
    H_part(i) = H_dis(q);
end
H_joint = H_dis(p);

MI = sum(H_part) - H_joint;

end