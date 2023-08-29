function [Aest,Best,Cest,Dest]= est_MaTS(Y,L,p,Z,Lz,pz,tol,maxit,printlevel)
% est_MaTS estimates the array pair (A,B) in a multi-component matrix valued AR(p)
% system. Estimation is performed via switching least squares to estimate
% the row and the column matrices. 
%
% SYNTAX:  [Aest,Best,Cest,Dest]= est_MaTS(Y,L,p)
%
% INPUT:   Y           ... T x M x N array of observations.
%          L           ... integer; lag length.
%          p           ... integer; number of terms per lag.
%          Z           ... T x Mz x Nz array of observations of exogenous vars.
%          Lz           ... integer; lag length for exo part.
%          pz           ... integer; number of terms per lag for exo part.
%          tol         ... real; tolerance for subsequent iterations
%          maxit       ... integer; maximum number of iterations 
%          printlevel ... integer; indicating, what information should be
%                          provided in the iterations. 
%
% OUTPUT:  Aest ... four dimensional array M x N x L x p. 
%          Best ... four dimensional array M x N x L x p.
%
% REMARKS: printlevel: + if >0: iteration number is printed.
%                       + if >1: info for each iteration is printed. 
% AUTHOR: dbauer, 6.6.2023. 
if nargin<8
    maxit = 100;
end
if nargin<7
    tol = 0.000000001;
end
if nargin<9
    printlevel = 0;
end
if nargin<6
    pz = 0;
end
if nargin<5
    Lz = 0;
end
if nargin<4
    Z = [];    
end

Leff = max(L,Lz); 

dims = size(Y);
T = dims(1)-Leff; % effective sample size
M = dims(2);
N = dims(3); 
 
if isempty(Z)
    exo = 0;
    Lz = 0;
    pz =0;
    Mz = 0;
    Nz = 0;
else
    exo = 1;
    dimsz = size(Z);
    Mz = dimsz(2);
    if length(dimsz)>2
        Nz = dimsz(3);
    else
        Nz = 1;
    end
end

% initialize randomly
Aest = randn(M,M,L,p)*0.1;
Best = randn(N,N,L,p)*0.1;
if exo
    Cest = randn(M,Mz,Lz,pz)*0.1;
    Dest = randn(N,Nz,Lz,pz)*0.1;
end 

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
    if exo 
        Cest_old = Cest;    
        Dest_old = Dest;
    end
    err_old = err; 
    % estimate (A,C) for given (B,D)
    % set up matrices
    for n=1:N
        % dependent variable
        vY(:,(n-1)*T+(1:T)) = squeeze(Y(Leff+(1:T),:,n))';
        % regressors
        for j=1:L
            for k=1:p
                for m=1:M
                    vX((j-1)*M*p+(k-1)*M+m,(n-1)*T+(1:T)) = squeeze(Best(n,:,j,k))*squeeze(Y(Leff-j+(1:T),m,:))';
                end
            end
        end

        if exo
            for j=1:Lz
                for k=1:pz
                    for m=1:Mz
                        vX(M*L*p+(j-1)*Mz*pz+(k-1)*Mz+m,(n-1)*T+(1:T)) = squeeze(Dest(n,:,j,k))*squeeze(Z(Leff-j+1+(1:T),m,:))';
                    end
                end
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
    if exo
        for j=1:Lz
            for k=1:pz
                Cest(:,:,j,k) = hatAC(:,M*L*p+(j-1)*Mz*pz+(k-1)*Mz+(1:Mz));
            end
        end
    end

    % estimate B for given A 
    % set up matrices
    for m=1:M
        % dependent variable
        vtY(:,(m-1)*T+(1:T)) = squeeze(Y(Leff+(1:T),m,:))';
        % regressors
        for j=1:L
            for k=1:p
                for n=1:N
                    vtX((j-1)*N*p+(k-1)*N+n,(m-1)*T+(1:T)) = squeeze(Aest(m,:,j,k))*squeeze(Y(Leff-j+(1:T),:,n))';
                end
            end
        end
        if exo
            for j=1:Lz
                for k=1:pz
                    for n=1:Nz
                        vtX(N*L*p+(j-1)*Nz*pz+(k-1)*Nz+n,(m-1)*T+(1:T)) = squeeze(Cest(m,:,j,k))*squeeze(Z(Leff-j+1+(1:T),:,n))';
                    end
                end
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
    if exo
        for j=1:Lz
            for k=1:pz
                Dest(:,:,j,k) = hatBD(:,N*L*p+(j-1)*Nz*pz+(k-1)*Nz+(1:Nz));
            end
        end
    end

    % renormalize 
    [Aest,Best] = norm_MaTS_syst(Aest,Best);
    if exo 
        [Cest,Dest] = norm_MaTS_syst(Cest,Dest);
    end
    % calculate difference in estimate 
    if exo 
        diff_exo = norm(Cest-Cest_old,"fro")+ norm(Dest-Dest_old,"fro");
    else 
        diff_exo= 0;
    end;
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
