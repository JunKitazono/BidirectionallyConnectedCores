function [ Cov_X, Cov_XY, Cov_Y ] = symAtoCovs( A, sigma_E )
%UNTITLED3 ���̊֐��̊T�v�������ɋL�q
%   �ڍא����������ɋL�q

N = size(A, 1);

Cov_X = sigma_E^2 * eye(N,N) / (eye(N,N)-A*A);
Cov_Y = Cov_X;
Cov_XY = sigma_E^2 * ( (eye(N,N)-A*A)\A );

end

