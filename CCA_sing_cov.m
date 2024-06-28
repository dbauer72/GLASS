function [th,A,K,C,D,Omega] = CCA_sing_cov(GammaT,n,q,kcol,krow,plots);
% CCA implements canonical correlations analysis subspace method for
% singular state space processes.
%
% SYNTAX: [th,a,b,c,d,k,Omega] = CCA_sing_cov(GammaT,n,kcol,krow,p);
%
% INPUTS:   GammaT ... NxNxkmax matrix of covariances up to lag kmax-1.
%           n ... integer; system order; if n=[] estimated using SVC.
%           q ... integer; dimension of singular input, if q<s. 
%           kcol, krow ... integers; past and future horizons for CCA.
%           plots ... indicator; if plots>0 singular values are plotted. 
%
% OUTPUTS:  th ... theta structure of estimated system.
%           (A,K,C) ... estimated state space system
%           Omega ... sxs estimated innovation variance matrix.
% 
% REMARK: CCA subspace algorithm (see Larimore 1987), no exogenous inputs; uses a regression framework
%
% dbauer, 30.1.2023

% first scale to improve numerical fit. 
SCALE = 1; %diag(1./sqrt(diag(squeeze(GammaT(:,:,1)))));

[N,~,kmax]=size(GammaT);

if kmax < krow+kcol
    disp('Not enough covariances. Reducing krow and kcol.')
    krow = floor(kmax/2);
    kcol = kmax - krow; 
end
for j=1:size(GammaT,3)
    GammaT(:,:,j) = SCALE*squeeze(GammaT(:,:,j))*SCALE;
end


% ---- three matrices -----
Wf2 = zeros(N*krow,N*krow);
for ja = 1:krow
    for jb = 1:krow 
        if (ja>= jb)
            Gj = squeeze(GammaT(:,:,ja-jb+1));
        else
            Gj = squeeze(GammaT(:,:,jb-ja+1))';
        end
        Wf2(N*(ja-1)+[1:N],N*(jb-1)+[1:N]) = Gj;
    end
end

Wp2 = zeros(N*kcol,N*kcol);
for ja = 1:kcol
    for jb = 1:kcol 
        if (ja>= jb)
            Gj = squeeze(GammaT(:,:,ja-jb+1))';
        else
            Gj = squeeze(GammaT(:,:,jb-ja+1));
        end
        Wp2(N*(ja-1)+[1:N],N*(jb-1)+[1:N]) = Gj;
    end
end

Hfp = zeros(N*krow,N*kcol);
for ja = 1:krow
    for jb = 1:kcol 
        Gj = squeeze(GammaT(:,:,ja+jb));
        Hfp(N*(ja-1)+[1:N],N*(jb-1)+[1:N]) = Gj;
    end
end

% ---- SVD ----
% do not use weight for future, as it might be singular. 
%Wf2 = (Yf*Yf')/T;
%Wf = inv(chol(Wf2)'); 
[u,s]=svd(Wf2);
ds = sqrt(diag(s));
ds(ds<0.0001)=0.0001;
is = diag(1./ds);
iWf =  u*is*u';

[u,s]=svd(Wp2);
ds = sqrt(diag(s));
ds(ds<0.0001)=0.0001;
s = diag(1./ds);

iWp =  u*s*u';

wbetaz = iWf * Hfp*iWp;
[U,S,V] = svd(wbetaz);

if plots   
    figure;
    plot(diag(S),'x');
    hold on;
end;

beta = V(:,1:n)'*iWp;
%beta = beta(:,1:n)\beta;


% ------ State variances -----
P = beta * Wp2 * beta'; 
iSCALE = inv(SCALE);

% ------ system matrices ------
Of = Hfp* beta'*inv(P); 
C = iSCALE*Of(1:N,:);

% ------ Omega -----

Omega = iSCALE*Wf2(1:N,1:N)*iSCALE - C*P*C';
[u,s,v] = svd(Omega); 
D= u(:,1:q); 

LamO = D'*Omega*D;

% ----- A and K -----
Gxe = [P,zeros(n,q);zeros(q,n),LamO];
beta1 = beta(:,1:N);
beta2 = [beta(:,(N+1):end),zeros(n,N)];

Gxp1xe = [(beta1 * Hfp(1:N,:) + beta2*Wp2)*beta',beta1*SCALE*D*LamO];
AK = Gxp1xe*inv(Gxe);

A = AK(:,1:n); 
K = AK(:,n+1:end);

%A = Of(1:(N*(krow-1)),:)\Of((N+1:N*krow),:);



% ----------   transformation to Echelon form
% generic neighbourhood
th = ss2ech_n(A,[K,zeros(n,N-q)],C);
th.D = D; 
th.B = th.K(:,1:q);

th.Omega = LamO; 

