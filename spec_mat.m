function Sig_est = spec_mat(y,freq,BT);
% estimation of smoothed periodogram 
%
% SYNTAX: Sig_est = spec_mat(y,freq,BT);
%
% INPUT: y ... real NxT matrix of data.
%        freq ... real Mx1 vector of frequencies for evaluation of spectral
%        estimate. Fractions of 2*pi. 
%        BT ... integer; smoothing. 
%
% OUTPUT:  Sig_est NxNxM array of estimated spectrum.
%
% REMARK: covariance based estimation, produces effectively a elementwise
% estimate.
%
% AUTHOR: dbauer, 27.1.2023. 

if nargin<3
    BT = length(freq);
end

[N,T]= size(y);

if BT > T/2
    BT = T/2;
    disp('Correct BT to T/2.');
end 
% detrend y
%for n=1:N
%    y(n,:)= detrend(y(n,:)',1)';
%end

% calculate covariances. 
Gammak = zeros(N,N,BT+1);
for k=1:BT+1
    Gammak(:,:,k) = y(:,1:(T-k+1))*y(:,k:(T))'/(2*pi*(T-k));
end

% calculate complex values.

F = length(freq);
Sig_est = zeros(N,N,F);
for f=1:F
    Sig_est(:,:,f)= squeeze(Gammak(:,:,1));
    om = exp(-freq(f)*2*pi*sqrt(-1));
    for k=1:BT
        Ga =squeeze(Gammak(:,:,k+1));
        Sig_est(:,:,f)= Sig_est(:,:,f)  + (1-(k/BT))*(Ga*om^k + Ga'*om^(-k));
    end
end
