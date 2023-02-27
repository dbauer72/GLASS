function [th,A,K,C,D,Omega] = CCA_sing(z,n,q,kcol,krow,plots);
% CCA implements canonical correlations analysis subspace method for
% singular state space processes.
%
% SYNTAX: [th,a,b,c,d,k,Omega] = CCA(z,n,kcol,krow,p);
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
Wf2 = Yf*Yf'/T; 
[u,s]=svd(Wf2);
ds = sqrt(diag(s));
ds(ds<0.000001)=0.000001;
is = diag(1./ds);
iWf =  u*is*u';


Hfp = Yf*Zp'/T; 
Wp2 = (Zp*Zp')/T;
[u,s]=svd(Wp2);
ds = sqrt(diag(s));
ds(ds<0.000001)=0.000001;
s = diag(1./ds);

iWp =  u*s*u';

wbetaz = iWf * Hfp*iWp;
[U,S,V] = svd(wbetaz);

% order estimation using SVC
if isempty(n)|(plots)
    s = diag(S);
    svc = s.^2 + log(T)/T*2*[1:length(s)]'*nz;
    nmax = length(s)-1;
    [minn,nhat] = min(svc(1:nmax+1));
    if (nhat>1)
       nhat = nhat-1;
    else
        nhat=1;
    end;
    
    nhat = max(nhat,1); % make sure that the estimated order is not too small.
end

if plots
    
    figure;
    plot(s,'x');
        hold on;
    plot([0,nmax],sqrt(1-log(T)/T)*[1,1],'m');
    chat = sum(s>sqrt(1-log(T)/T));
    title(sprintf('Singular values: c: %d, SVC: %d',chat,nhat));
    plot(nhat,S(nhat,nhat),'ro');
    for c=1:chat,
        plot(c,S(c,c),'mo');
    end;
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
D = D(:,1:q);
eps = D'*res;
% ------ residuals will be singular -----

AK = x(:,2:end)/[x(:,1:end-1);eps(:,1:end-1)];
A = AK(:,1:n);
K = AK(:,n+1:end);


% ----------   transformation to Echelon form
% generic neighbourhood
th = ss2ech_n(A,[K,zeros(n,nz-q)],C);
th.D = D; 
th.B = th.K(:,1:q);
th.Omega = LamO(1:q,1:q); 

