function [th,A,B,C,D,K,Omega] = CCAX(z,s,n,kcol,krow,plots,D0);
% CCAX implements canonical correlations analysis subspace method including exogenous inputs.
%
% SYNTAX: [th,A,B,C,D,K,Omega] = CCAX(z,s,n,kcol,krow,plots,D0);
%
% INPUTS:   z ... Tx (s+m) matrix of observations.
%           n ... integer; system order; if n=[] estimated using SVC.
%           kcol, krow ... integers; past and future horizons for CCA.
%           plots ... indicator; if plots>0 singular values are plotted. 
%           D0   ... indicator; if D0>0 (default: 1), then D is estimated, else it is
%                             seen as zero. 
%
% OUTPUTS:  th ... theta structure of estimated system.
%           (A,B,C,D,K) ... estimated state space system
%           Omega ... sxs estimated innovation variance matrix.
% 
% REMARK: CCA subspace algorithm (see Larimore 1987), with exogenous
% inputs; uses a regression framework.
%
% dbauer, 26.4.2023
if nargin<7
    D0=0;
end

if D0 ~= 0
    D0 = 1;
end

[T,nz]=size(z);
m = nz-s; 
if (m<0)
    error('Wrong dimension, not space left for exogenous inputs!');
end;
% ------ data Hankel matrices -----
while krow*nz> (T-krow-kcol)
    krow = krow-1;
end

for i=1:krow
   Yf((i-1)*s+[1:s],:) = z(kcol+i+[0:T-krow-kcol],1:s)';
   Uf((i-1)*m+[1:m],:) = z(kcol+i+[0:T-krow-kcol],(s+1):nz)';
end;
for i=1:kcol
   Zp((kcol-i)*nz+[1:nz],:) = z(i+[0:T-krow-kcol],1:nz)';
end;
 
% regress out Uf -> this can be done iteratively, only regressing out the
% relevant terms. 
Yfperp = Yf;
for i=(1+D0):krow
    y = Yf((i-1)*s+[1:s],:)';
    u = Uf(1:((i-D0)*m),:)';
    y = y - u* (u\y);
    Yfperp((i-1)*s+[1:s],:) = y';
end

Zpperp = Zp;

% ------ regressions to obtain hat beta ----
betaz = Yfperp/Zpperp;

% ---- SVD ----
Wf2 = (Yfperp*Yfperp')/T;
Wf = inv(chol(Wf2)');

Wp2 = (Zpperp*Zpperp')/T;
if min(eig(Wp2))<10^(-6)
    Wp2 = Wp2 + eye(size(Wp2))*10^(-6);
end
Wp =  chol(Wp2)';

[U,S,V] = svd(Wf*betaz*Wp);

% order estimation using SVC
if isempty(n)|(plots)
    dS = diag(S);
    svc = dS.^2 + log(T)/T*2*[1:length(dS)]'*nz;
    nmax = length(dS)-1;
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
    plot(dS,'x');
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
beta = V(:,1:n)'*inv(Wp);
beta = beta(:,1:n)\beta;

% ------ system matrices ------
x = beta*Zp;
if (D0==0)
    CD = Yf(1:s,:)/[x;Uf(1:m,:)];
    C = CD(:,1:n);
    D = CD(:,n+1:end);
else 
    CD = Yf(1:s,:)/[x];
    C = CD(:,1:n);
    D = zeros(s,m);
end

res = Yf(1:s,:) - C*x-D*Uf(1:m,:);
ABK = x(:,2:end)/[x(:,1:end-1);Uf(1:m,1:end-1);res(:,1:end-1)];
A = ABK(:,1:n);
B = ABK(:,n+(1:m));
K = ABK(:,(n+m+1):end);
Omega = (res*res')/T;

% ----------   transformation to Echelon form
% generic neighbourhood
th = ss2ech_n(A,[B,K],C);
th.D = D;
th.Omega = Omega; 
