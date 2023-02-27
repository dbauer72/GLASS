function Ce = condense_square_root(LTs,plots);
% condense square root takes the matrix square roots in the array LTs and
% extracts a rank r decomposition of the summed square matrix. 
%
% SYNTAX: Ce = condense_square_root(LTs);
%
% INPUT:  LTs   ... NxqxM real matrix; reduced rank matrix square roots
%         plots ... indicator for plotting the singular values. 
%
% OUTPUT: Ce  Nxq real matrix; dominant directions.
%
% AUTHOR: dbauer, 27.1.2023.
if (nargin<2)
    plots=0;
end,

[N,r,M] = size(LTs);

% aggregate 
Sig = zeros(N,N); 
for m=1:M
    L = squeeze(LTs(:,:,m));
    Sig = Sig + L*L';
end

% decompose. 
[U,S] = svd(Sig);

Ce = U(:,1:r);
if plots
    plot(log(diag(S)),'x');
end