function R = cal_cov(th,k);
% calculate the covariance sequence corresponding to the system th.
%
% assumes that the system is stable. 
%
% author: dbauer, 19.12.2024. 

A= th.A;
C = th.C;
K = th.K; 
Om = th.Omega; 

[s,n]=size(C);

% calculate P 
vQ = K*th.Omega*K';
P = reshape( inv(eye(n^2) - kron(A,A))*vQ(:),n,n);
M = A*P*C' + K*Om; 

% define sequence 
R = zeros(s,s,k+1); 
R(:,:,1) = C*P*C'+Om; 

for j=1:k
    R(:,:,j+1) = C*A^(j-1)*M;
end;
