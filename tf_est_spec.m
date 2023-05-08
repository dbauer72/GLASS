function [the,GammaT] = tf_est_spec(LTs,freq,q,n,k)
% tf_est_spec provides direct estimation of a state space system for 
% low rank spectral factors LTs evaluated at frequencies freq.
% The state space system has order n and is estimated using CCA_sing using
% a maximum of k lags. 
%
% SYNTAX: the = tf_est_spec(LTs_chi,freq,n,k)
%
% INPUT:    LTs  ... NxrxF real array; spectral factors at frequencies.
%           freq ... F real vector; frequencies at which spectrum is evaluated.
%           n    ... integer; order of state space system
%           k    ... integer; max lag length for subspace algorithm.
%
% OUTPUT:  the ... theta structure. 
%
% AUTHOR: dbauer, 30.1.2023.

[N,r,F] = size(LTs);

if length(freq) ~= size(LTs,3)
    disp('Number of frequencies does not match size of LTs!')
    return;
end

% extract covariance matrix
GammaT = zeros(N,N,2*k); % subspace algorithm needs 2*k-1 covariances. 
om = exp(sqrt(-1)*freq*2*pi);
df = 2*pi*(freq(2)-freq(1));
for f=1:F
    Lt  = squeeze(LTs(:,:,f));
    Omegaf = (Lt*Lt');
    for jj = 1:(2*k)
        %% 
        GammaT(:,:,jj)=GammaT(:,:,jj) + real(Omegaf*om(f)^(1-jj) + Omegaf.'*om(f)^(jj-1));
    end
end

GammaT = GammaT*df; 
% use CVA to estimate the system
[the] = CCA_sing_cov(GammaT,n,q,k,k,1);

