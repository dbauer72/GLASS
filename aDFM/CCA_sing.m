function [th,A,K,C,D,Omega,nhat,nhat2,se,x] = CCA_sing(z,n,q,kcol,krow,plots);
% CCA implements canonical correlations analysis subspace method for
% singular state space processes.
%
% SYNTAX: [th,A,K,C,D,Omega] = CCA_sing(z,n,q,kcol,krow,plots);
%
% INPUTS:   z ... Tx s matrix of observations.
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
% dbauer, 27.10.2019


[T,nz]=size(z);

% ------ weigh matrix to get rid of integration annulations 
y = z([2:T],1:nz);
dy = y-z([1:(T-1)],1:nz);

WW = dy'*dy/(T-1);
CW = eye(nz); %chol(WW);
iCW = inv(CW);

z(:,1:nz) = z(:,1:nz)*iCW;
% ------ data Hankel matrices -----
while krow*nz> (T-krow-kcol)
    krow = krow-1;
end

for i=1:krow
   Yf((i-1)*nz+[1:nz],:) = z(kcol+i+[0:T-krow-kcol],1:nz)';
end;
for i=1:kcol
   Zp((kcol-i)*nz+[1:nz],:) = z(i+[0:T-krow-kcol],1:nz)';
end;
      
% ------ regressions to obtain hat beta ----



% ---- SVD ----
% do not use weight for future, as it might be singular. 
%Wf2 = (Yf*Yf')/T;
%Wf = inv(chol(Wf2)');
tol = 0.000001;
Wf2 = Yf*Yf'/T; 
[uf,sf]=svd(Wf2);
ds = sqrt(diag(sf));
ds(ds<tol)=tol;
is = diag(1./ds);
iWf_cca =  uf*is*uf';
iWf = eye(size(Yf,1));

Hfp = Yf*Zp'/T; 
Wp2 = (Zp*Zp')/T;
[up,sp]=svd(Wp2);
ds = sqrt(diag(sp));
%ds = (diag(sp));

ds(ds<tol)=tol;
se = diag(1./ds);

iWp =  up*se*up';

wbetaz = iWf * Hfp*iWp;
[U,S,V] = svd(wbetaz);

% order estimation 
% %using SVC
%if isempty(n)|(plots)
    s = diag(S);
    svc = s.^2 + log(T)/T*2*[1:length(s)]'*nz;
    nmax = length(s)-1;
    [minn,nhat] = min(svc(1:nmax+1));
    if (nhat>1)
       nhat = nhat-1;
    else
        nhat=1;
    end;
    
    nhat2 = max(nhat,1); % make sure that the estimated order is not too small.
%end

% alternative: use can. corr. 
betaz = iWf_cca * Hfp*iWp;
[Ucca,Scca,Vcca] = svd(betaz);
se = diag(Scca); 
nhat2 = sum(se>sqrt(1-log(T)/T));
    
if plots
    
    figure;
    plot(se,'x');
        hold on;
    plot([0,nmax],sqrt(1-log(T)/T)*[1,1],'m');
    title(sprintf('Singular values: c: %d, SVC: %d',nhat,nhat2));
    if nhat>0
        plot(nhat,Scca(nhat,nhat),'ro');
        for c=1:nhat,
            plot(c,Scca(c,c),'mo');
        end;
    end
end;

if isempty(n)
    n = nhat;
end;
beta = V(:,1:n)'*iWp;
%beta = beta(:,1:n)\beta;

% ------ system matrices ------
x = beta*Zp;

C = Yf(1:nz,:)/x;
res = Yf(1:nz,:) - C*x;
Omega = (res*res')/T;

[D,LamO]= svd(Omega); 
if isempty(q)
    q = select_q(Omega,T);
end

D = D(:,1:q);
eps = D'*res;
% ------ residuals will be singular -----
AK = x(:,2:end)/[x(:,1:end-1);res(:,1:end-1)];
A = AK(:,1:n);
K = AK(:,n+1:end)*D;


% ----------   transformation to Echelon form
% generic neighbourhood
th = ss2ech_n(A,[K,zeros(n,nz-q)],CW'*C);
th.D = CW'*D; 
th.B = th.K(:,1:q);
th.Omega = LamO(1:q,1:q); 

