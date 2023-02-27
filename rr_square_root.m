function [LTs,sigmas] = rr_square_root(Sig_est,q,plots);
% provides the square roots of the set of matrices in the dominant q
% dimensions.
%
% SYNTAX: LTs = rr_square_root(Sig_est,q,plots);
%
% INPUT: Sig_est ... NxNxM real array of spectrum estimates at different
% frequencies. 
%         q ... integer; number of components. 
%         plots ... indicator, if aggregated singular values should be
%         plotted. 
%
% OUTPUT: LTs ... NxqxM real matrix of dominant directions. 
%         sigmas ... NxM real matrix of estimated singular values. 
%
% AUTHOR: dbauer, 27.1.2023. 
if nargin<3
    plots = 0;
end
[N,~,M] = size(Sig_est);

LTs = zeros(N,q,M);
sigmas = zeros(N,M); 

for m=1:M
    sige = squeeze(Sig_est(:,:,m));
    [u,s] = svd(sige); 
    sigmas(:,m) = sqrt(diag(s));
    LTs(:,:,m)=u(:,1:q)*diag(sigmas(1:q,m)); 
end

if plots
    figure
    plot(sigmas');
end
