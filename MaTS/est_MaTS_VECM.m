function [Aest,Best,alphaest,betaest,Cest,Dest]= est_MaTS_VECM(Y,r,L,p,tol,maxit,printlevel)
% est_MaTS_VECM estimates the array pair (A,B) in a multi-component matrix valued AR(p)
% system as well as alpha beta' in the vecotized VECM equation. 
% Estimation is performed via switching least squares to estimate
% the row and the column matrices. 
%
% SYNTAX:  [Aest,Best,alphaest,betaest]= est_MaTS_VECM(Y,r,L,p)
%
% INPUT:   Y           ... T x M x N array of observations.
%          r           ... pair of integers; cointegrating rank and number of terms.  
%          L           ... integer; lag length.
%          p           ... integer; number of terms per lag.
%          tol         ... real; tolerance for subsequent iterations
%          maxit       ... integer; maximum number of iterations 
%          printlevel ... integer; indicating, what information should be
%                          provided in the iterations. 
%
% OUTPUT:  Aest ... four dimensional array M x N x L x p. 
%          Best ... four dimensional array M x N x L x p.
%          alphaest ... real M*N x r(1) matrix; estimated loading matrix.
%          betaest  ... real M*N x r(1) matrix; cointegrating relations. 
%
% REMARKS: printlevel: + if >0: iteration number is printed.
%                       + if >1: info for each iteration is printed. 
% AUTHOR: dbauer, 27.6.2024. 
if nargin<6
    maxit = 100;
end
if nargin<5
    tol = 0.000000001;
end
if nargin<7
    printlevel = 0;
end

dims = size(Y);
T = dims(1)-L-1; % effective sample size
M = dims(2);
N = dims(3); 
 
dY = Y(2:end,:,:) - Y(1:end-1,:,:);

% initialize randomly
Aest = randn(M,M,L,p)*0.1;
Best = randn(N,N,L,p)*0.1;

% define Pi vY_{t-1} as exogenous term 
rr = r(1);
Jr = r(2);
Lz = 1;
pz = Jr; 
Mz = M;
Nz = N;

% for initial values no cointegrating rank is used.
Cest = randn(M,M,1,Jr)*0.1;
Dest = randn(N,N,1,Jr)*0.1;

Ym1 = Y(1:end-1,:,:);


% predefined for numerical speed. 
% Warning: These matrices might get really large!
vY = zeros(M,N*T);
vtY = zeros(N, M*T);
vX= zeros(M*L*p+Mz*Lz*pz,N*T);
vtX = zeros(N*L*p+ Nz*Lz*pz,M*T);

% iteration 
err = norm(Y,"fro");
it = 0;
change = 1;
while ((change>tol)&&(it<maxit))
    % update 
    Aest_old = Aest;
    Best_old = Best; 

    Cest_old = Cest;    
    Dest_old = Dest;

    err_old = err; 
    % estimate (A,C) for given (B,D)
    % set up matrices
    for n=1:N
        % dependent variable
        vY(:,(n-1)*T+(1:T)) = squeeze(dY(L+(1:T),:,n))';
        % regressors
        for j=1:L
            for k=1:p
                for m=1:M
                    vX((j-1)*M*p+(k-1)*M+m,(n-1)*T+(1:T)) = squeeze(Best(n,:,j,k))*squeeze(dY(L-j+(1:T),m,:))';
                end
            end
        end

        for k=1:pz
            for m=1:Mz
                vX(M*L*p+(k-1)*Mz+m,(n-1)*T+(1:T)) = squeeze(Dest(n,:,1,k))*squeeze(Ym1(L+(1:T),m,:))';
            end
        end
    end
    % calculate LS fit 
    hatAC = vY / vX;
    % fill in estimates
    for j=1:L
        for k=1:p
            Aest(:,:,j,k) = hatAC(:,(j-1)*M*p+(k-1)*M+(1:M));
        end
    end


    for k=1:pz
        Cest(:,:,1,k) = hatAC(:,M*L*p+(k-1)*Mz+(1:Mz));
    end

    % estimate B for given A 
    % set up matrices
    for m=1:M
        % dependent variable
        vtY(:,(m-1)*T+(1:T)) = squeeze(dY(L+(1:T),m,:))';
        % regressors
        for j=1:L
            for k=1:p
                for n=1:N
                    vtX((j-1)*N*p+(k-1)*N+n,(m-1)*T+(1:T)) = squeeze(Aest(m,:,j,k))*squeeze(dY(L-j+(1:T),:,n))';
                end
            end
        end

        for k=1:pz
            for n=1:Nz
               vtX(N*L*p+(k-1)*Nz+n,(m-1)*T+(1:T)) = squeeze(Cest(m,:,1,k))*squeeze(Ym1(L+(1:T),:,n))';
            end
        end
    end
    % calculate LS fit 
    hatBD = vtY / vtX;
    % fill in estimates
    for j=1:L
        for k=1:p
            Best(:,:,j,k) = hatBD(:,(j-1)*N*p+(k-1)*N+(1:N));
        end
    end

    for j=1:Lz
        for k=1:pz
           Dest(:,:,1,k) = hatBD(:,N*L*p+(k-1)*Nz+(1:Nz));
        end
    end


    % renormalize 
    [Aest,Best] = norm_MaTS_syst(Aest,Best);

    [Cest,Dest] = norm_MaTS_syst(Cest,Dest);

    % calculate difference in estimate 
 
    diff_exo = norm(Cest-Cest_old,"fro")+ norm(Dest-Dest_old,"fro");
    diff_syst = norm(Aest-Aest_old,"fro") + norm(Best-Best_old,"fro") + diff_exo;
    err = norm(vtY-hatBD*vtX,"fro");
    it = it+1;

    change = max((err-err_old),diff_syst);
    % print out results, if wanted 
    if printlevel>0
        fprintf("Iteration #%d (of max %d) \n",it,maxit);
        if printlevel>1
            fprintf("Function value: %f. Improvement: %f\n",err^2/(T*M*N), (err_old-err)/err_old);
            fprintf("Step size: %f \n\n\n", diff_syst);
        end
    end
    % done -> iterate 
end

% estimate cointegrating relation using SVD of Pi. 

Pi = AB_to_vAj(Cest,Dest);
[u,s,v] = svd(Pi);

alphaest = u(:,1:rr)*s(1:rr,1:rr);
betaest = v(:,1:rr);

% normalize 
[alphaest,betaest] = norm_MaTS_VECM(alphaest,betaest,M,N,r);

