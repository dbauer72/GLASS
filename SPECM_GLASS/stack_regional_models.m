function [th_full,ve_full,x_full] = stack_regional_models(ths,indices,Wstar,Ww,y_full,Pbull);
% This functions stacks the regional models and compiles a joint model to
% evaluate the GLASS approach.
%
% SYNTAX: [th_stacked,ve_full] = stack_regional_models(ths,indices,Wstar,Ww,y_full);
%
% INPUT:
%           ths ... vector of regional models in innovation form.
%           indices... Nx3 matrix of integers, row i contains the index for the i-th regional model.  
%           Wstar ... N vector of weighting matrices for star variables
%           Ww   ...  d x N_full (total sum of vars)
%           y_full ... N_full x T real matrix of observations. 
%           P_bull ... integer indicating whether the PE approach should be
%                       used.  
%
% OUTPUT: 
%           th_stacked ... theta model for stacked process.
%           ve_full    ... N_full x T real matrix of estimated innovations.
%           x_full     ... n_f x T real matrix of estimated states. 
%
% AUTHOR: dbauer, 5.8.2024. 

[N_full,T]=size(y_full);
N = size(indices,1);
csi = [0,;cumsum(indices(:,1))];

N_full = csi(end);
th_stacked = theta_urs();
ve_full = y_full*0; 


nis = zeros(N,1);
for jn=1:N
    nis(jn)= size(ths(jn).A,1);
end

cnis = [0;cumsum(nis)];
x_full = zeros(cnis(end),T); 

wt = Ww*y_full; % dominant variables. 

% initialize matrices 
n_full = cnis(end);
A_full = zeros(n_full,n_full);
C_full = zeros(N_full,n_full);
K_full = zeros(n_full,N_full);

D_full = zeros(N_full,N_full); 

for jn = 1:N % cycle over regions
    thi = ths(jn);
    ni = size(thi.A,1);
    yi = [y_full(csi(jn)+1:csi(jn+1),:)];
    yistar = Wstar{jn}*y_full;

    sist = size(yistar,1);
    si = indices(jn,1);
    yi_full = [yi;yistar;wt];
    % provide results 
    parc_star = th2param_RM(thi,indices(jn,:),1);
    paromi = extr_lowtri(thi.Omega);
    param_star = [paromi(:);parc_star(:)];
    [ql,tresi] = cal_quasi_like_RM(param_star,yi_full',ni,indices(jn,:),Pbull);

    % orthogonalize epsilons. 
    Omegai = tresi'*tresi/T; 
    Di = Omegai(1:si,(si+1):end)*inv(Omegai((si+1):end,(si+1):end));
    tildeepsi = tresi(:,1:si)-tresi(:,(si+1):end)*Di';

    % fill in residuals 
    ve_full(csi(jn)+1:csi(jn+1),:) = tildeepsi'; 
    % filter to obtain estimate of the state. 
    x_full(cnis(jn)+1:cnis(jn+1),:) = ltitr(thi.A-thi.K*thi.C,thi.K,yi_full',zeros(ni,1))';

    % calculate converted system: 
    tCi = thi.C(1:si,:) - Di*thi.C((si+1):end,:);
    tKistar = thi.K(:,1:si)*Di + thi.K(:,(si+1):end);
    tAi = thi.A - tKistar*thi.C((si+1):end,:) - thi.K(:,1:si)*tCi;
    tKi = [thi.K(:,1:si),tKistar];

    % fill into the final matrices 
    Wi = [Wstar{jn};Ww];
    D_full(csi(jn)+1:csi(jn+1),csi(jn)+1:csi(jn+1)) = eye(si);
    D_full(csi(jn)+1:csi(jn+1),:) = D_full(csi(jn)+1:csi(jn+1),:)-Di*Wi;

    A_full(cnis(jn)+1:cnis(jn+1),cnis(jn)+1:cnis(jn+1))= tAi;
    C_full(csi(jn)+1:csi(jn+1),cnis(jn)+1:cnis(jn+1)) = tCi; 
    K_full(cnis(jn)+1:cnis(jn+1),csi(jn)+1:csi(jn+1)) = thi.K(:,1:si);
    K_full(cnis(jn)+1:cnis(jn+1),:) = K_full(cnis(jn)+1:cnis(jn+1),:)+tKistar*Wi;    
end

% inverse system: D_full y_full = C_full x_full + tres; 
%     x_full(t+1) = A_full x_full(t) + K_full y(t). 

sev = svd(D_full); 
if (min(sev)<0.000001)
    disp('D_full is close to singular!')
end

tildeC_full = -inv(D_full)*C_full; 

% invert system. 
th_full.A = A_full - K_full*tildeC_full; 
th_full.K = K_full;
th_full.C = -tildeC_full; 
th_full.B = zeros(n_full,0);

Omega_full = ve_full*ve_full'/T;
th_full.Omega = Omega_full; 




