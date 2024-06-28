function [thetaest,Aest,Best,Vest,sdA,sdB,sPsi,par_init,Ai,Bi]= qMLE_MaTS(Y,L,p)
% qMLE_MaTS estimates the array pair (A,B) in a multi-component matrix valued p term AR(L)
% system. Estimation is performed via quasi maximum likelihood maximization using the Gaussian likelihood.
%
%
% SYNTAX:  [Aest,Best,Cest,Dest]= est_MaTS(Y,L,p)
%
% INPUT:   Y           ... T x M x N array of observations.
%          L           ... integer; lag length.
%          p           ... integer; number of terms per lag.
%
% OUTPUT:  thetaest ... real d vector; estimated parameter vector
%          Aest     ... four dimensional array M x N x L x p. 
%          Best     ... four dimensional array M x N x L x p.
%          Vest     ... real d x d matrix; estimated variance matrix for thetaest.
%
% REMARKS: identification for each lag is achieved using the QR
% decomposition of the Phi transformed sum of of
% vectorized terms. 
%
% AUTHOR: dbauer, 22.8.2023. 

%%%          Z           ... T x Mz x Nz array of observations of exogenous vars.
%%%          Lz           ... integer; lag length for exo part.
%%%          pz           ... integer; number of terms per lag for exo part.
%%%          tol         ... real; tolerance for subsequent iterations
%%%          maxit       ... integer; maximum number of iterations 
%%%          printlevel ... integer; indicating, what information should be
%%%                          provided in the iterations. 


% get dimensions 
[T,M,N] = size(Y);

% get initial estimate
[Ai,Bi]= est_MaTS(Y,L,p);

% calucate corresponding parameter vector
par_init = MaTS2param(Ai,Bi);

% non-linear optimization
options = optimoptions('fminunc','display','iter','MaxIterations',1000,'Algorithm','trust-region','SpecifyObjectiveGradient',true);
options.MaxFunctionEvaluations = 1000;

[thetaest,exitflag] = fminunc(@(x) cal_crit_MaTS(Y,x,L,p),par_init,options);

% prepare output 
[crit,grad,Hess,JV] = cal_crit_MaTS(Y,thetaest,L,p);
iH = inv(Hess); 
Vest = iH * JV * iH;

% calculate AR coefficient matrices 
[Aest,Best,dA,dB] = param2MaTS(thetaest,M,N,L,p);

% calculate standard deviations for (A,B) coefficient matrices
sdA = Aest*0;
sdB = Best*0;

for mr=1:M
    for mc=1:M
        for l=1:L
            for jr=1:p
                da = squeeze(dA(mr,mc,l,jr,:));
                sdA(mr,mc,l,jr)= sqrt(da' * Vest* da);
            end
        end
    end
end
for nr=1:N
    for nc=1:N
        for l=1:L
            for jr=1:p
                db = squeeze(dB(nr,nc,l,jr,:));
                sdB(nr,nc,l,jr)= sqrt(db' * Vest* db);
            end
        end
    end
end

% calculate standard deviations for AR matrices 
dvPsi = zeros(M*N,M*N,L,length(thetaest));
for d=1:length(thetaest)
    for l=1:L
        for jr=1:p
            dvPsi(:,:,l,d) = dvPsi(:,:,l,d) + kron(squeeze(Best(:,:,l,jr)),squeeze(dA(:,:,l,jr,d))) + kron(squeeze(dB(:,:,l,jr,d)),squeeze(Aest(:,:,l,jr)));
        end
    end
end

sPsi = zeros(N*M,N*M,L); % standard deviations for all entries in all lags. 
for row=1:N*M
    for col=1:N*M
        for l=1:L
            dPsi_row = squeeze(dvPsi(row,col,l,:));
            sPsi(row,col,l)= sqrt(dPsi_row' * Vest* dPsi_row);
        end
    end
end


