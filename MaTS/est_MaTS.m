function [Aest,Best]= est_MaTS(Y,L,p,tol,maxit,printlevel)
% est_MaTS estimates the array pair (A,B) in a multi-component matrix valued AR(p)
% system. Estimation is performed via switching least squares to estimate
% the row and the column matrices. 
%
% SYNTAX:  [Aest,Best]= est_MaTS(Y,L,p)
%
% INPUT:   Y           ... T x M x N array of observations.
%          L           ... integer; lag length.
%          p           ... integer; number of terms per lag.
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
if nargin<5
    maxit = 100;
end
if nargin<4
    tol = 0.000000001;
end
if nargin<5
    printlevel = 0;
end
dims = size(Y);
T = dims(1)-L; % effective sample size
M = dims(2);
N = dims(3); 

% initialize randomly
Aest = randn(M,M,L,p)*0.1;
Best = randn(N,N,L,p)*0.1;

% predefined for numerical speed. 
% Warning: These matrices might get really large!
vY = zeros(M,N*T);
vtY = zeros(N, M*T);
vX= zeros(M*L*p,N*T);
vtX = zeros(N*L*p,M*T);


% iteration 
err = norm(Y,"fro");
it = 0;
change = 1;
while ((change>tol)&&(it<maxit))
    % update 
    Aest_old = Aest;
    Best_old = Best; 
    err_old = err; 
    % estimate A for given B
    % set up matrices
    for n=1:N
        % dependent variable
        vY(:,(n-1)*T+(1:T)) = squeeze(Y(L+(1:T),:,n))';
        % regressors
        for j=1:L
            for k=1:p
                for m=1:M
                    vX((j-1)*M*p+(k-1)*M+m,(n-1)*T+(1:T)) = squeeze(Best(n,:,j,k))*squeeze(Y(L-j+(1:T),m,:))';
                end
            end
        end
    end
    % calculate LS fit 
    hatA = vY / vX;
    % fill in estimates
    for j=1:L
        for k=1:p
            Aest(:,:,j,k) = hatA(:,(j-1)*M*p+(k-1)*M+(1:M));
        end
    end

    % estimate B for given A 
    % set up matrices
    for m=1:M
        % dependent variable
        vtY(:,(m-1)*T+(1:T)) = squeeze(Y(L+(1:T),m,:))';
        % regressors
        for j=1:L
            for k=1:p
                for n=1:N
                    vtX((j-1)*N*p+(k-1)*N+n,(m-1)*T+(1:T)) = squeeze(Aest(m,:,j,k))*squeeze(Y(L-j+(1:T),:,n))';
                end
            end
        end
    end
    % calculate LS fit 
    hatB = vtY / vtX;
    % fill in estimates
    for j=1:L
        for k=1:p
            Best(:,:,j,k) = hatB(:,(j-1)*N*p+(k-1)*N+(1:N));
        end
    end

    % renormalize 
    [Aest,Best] = norm_MaTS_syst(Aest,Best);
    
    % calculate difference in estimate 
    diff_syst = norm(Aest-Aest_old,"fro") + norm(Best-Best_old,"fro");
    err = norm(vtY-hatB*vtX,"fro");
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
