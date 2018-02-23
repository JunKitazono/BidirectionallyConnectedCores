function [p_past,joint,p_present] = est_prior_joint(X,N_s,tau)

%-------------------------------------------------------------------------------------------------
% PURPOSE: estimate probability distributions of past (p_past), present (p_present) and 
% joint distribution of past and present (joint) from discretized time series data X 
%
% INPUTS:
%   X: discretized time series data in form units x time
%   N_s: the number of states in each unit
%   tau: time lag between past and present
%
% OUTPUTS:
%   p_past: probability distribution of past state (X^t-tau)
%   joint: joint distribution of X^t (present) and X^(t-\tau) (past)
%   p_present: probability distribution of present state (X^t)
%
%-------------------------------------------------------------------------------------------------
%
% Masafumi Oizumi, 2018

%%

N = size(X,1);
T = size(X,2);
p_past = zeros(N_s^N,1);
p_present = zeros(N_s^N,1);
joint = zeros(N_s^N,N_s^N);

for t=1: T-tau
    past_s = baseMto10(X(:,t) - 1, N_s) + 1; % past state
    p_past(past_s) = p_past(past_s) + 1;
    present_s = baseMto10(X(:,t+tau) - 1, N_s) + 1; % present state
    p_present(present_s) = p_present(present_s) + 1;
    joint(present_s,past_s) = joint(present_s,past_s) + 1;
end

p_past = p_past/(T-tau);
p_present = p_present/(T-tau);
joint = joint/(T-tau);

end