function [crit,ve2] = cal_crit_MaTS_VECM(Y,theta,r,L,p)
% cal_crit_MaTS_VECM calculates the concentrated (w.r.t. the parameters for the variance) 
% Gaussian likelihood criterion function for the VECM scenario.
%
% SYNTAX: [crit,ve] = cal_crit_MaTS(Y,theta,M,N,L,p)
%
% INPUTS: Y       ... T x M x N real matrix of observations.
%         theta   ... d vector of parameters
%         M       ... integer; 
%         N       ... integer
%         r       ... pair of integers; cointegrating rank and number of terms.  
%         L       ... integer; number of lags
%         p       ... integer; number of terms per lag.
%
% OUTPUTS: crit   ... real, value of criterion function.
%          ve     ... T x M*N real matrix of vectorized residuals.
%
% REMARK: no derivatives calculated as of now. 
%
% AUTHOR: dbauer, 27.6.2024.

% get dimensions.
[T,M,N] = size(Y); 
nth = length(theta);
grad= zeros(nth,1);
Hess = zeros(nth,nth); 

% calculate AR matrices from parameters 
[GamA,GamB,alpha,beta] = param2MaTS_VECM(theta,M,N,r,L,p);
[vec_poly] = vectorized_syst(GamA,GamB);


% calculate residuals 
ve = zeros(T,M*N); 
% initial vallues for residuals set to NaN
ve(1:L,:) = NaN;  

% alternative: vectorize
vY = reshape(Y,T,M*N);
dvY = [vY(1,:);vY(2:end,:)-vY(1:(end-1),:)];
vYm1 = [NaN(1,M*N);vY(1:(end-1),:)];


ve2 = dvY((L+1):T,:)-vYm1((L+1):T,:)*beta*alpha';

for l=1:L
    ve2 = ve2 - dvY((L+1-l):(T-l),:)*squeeze(vec_poly(:,:,l))';
end

ve2(isnan(ve2))= 10^6;



Omegahat = ve2(L+1:end,:)'*ve2(L+1:end,:)/(T-L);

crit = log(det(Omegahat)+.001);

% if (any(abs(Omegahat(:))>10^6))
%     crit = log(max(Omegahat(:)));
% else 
%     % calculate criterion value 
%     crit = log(det(Omegahat));
% end







