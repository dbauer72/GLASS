function [Ae,Be,Ce,De,Omegae,Lambdae,Pse] = est_EM(y,xs,Ps,Pst,r); 
%
%

[s,T]=size(y); 
[n,~]= size(xs);


% cumulate Ps. 
Gam0 = sum(Ps(:,:,1:T),3);
Gam1 = sum(Pst(:,:,1:T),3);
Gamp = sum(Ps(:,:,2:(T)),3);

% regression matrices 
XtX = (xs(:,1:T)*xs(:,1:T)'+Gam0)/T;
Xty = (xs(:,1:T)*y(:,1:T)')/T; 
yty = (y(:,1:T)*y(:,1:T)'); 

% ZtZ, ZtX
ZtZ  = [yty, y(:,1:T)*xs(:,2:end)';xs(:,2:end)*y(:,1:T)',xs(:,2:end)*xs(:,2:end)'+Gamp]/T;
ZtX = [y(:,1:T)*xs(:,1:T)';xs(:,2:end)*xs(:,1:T)'+Gam1]/T;

% for C 
CAi = [y*xs(:,1:T)';xs(:,2:end)*xs(:,1:T)'+Gam1]*inv(XtX)/T;

LCi = CAi(1:s,:);
Ai = CAi((s+1):end,:);

% check stability and regularize, if needed.
[ev] = eig(Ai);
dev = abs(ev); 
mdev = max(dev); 
if (mdev>.99)
    Ai = Ai/mdev*0.99;
end

SigOmega = ZtZ- CAi*XtX*CAi';

% initial estimate for B,D, Psi. 
[U,S,V]=svd(SigOmega);
BDi = U(:,1:r)*sqrt(S(1:r,1:r));
[Q,R] = qr(BDi(1:r,1:r)');
BDi = BDi*Q; 
LDi = BDi(1:s,:);
Bi = BDi((s+1):end,:);


% extract estimator for Lambda.
CDi = [LCi,LDi];
[Q,R]=qr(CDi);

Lambdae = Q(:,1:r)*sqrt(s);
Lambdadagger = Q(:,1:r)'/sqrt(s);

Ci = Lambdadagger*LCi;
Di = Lambdadagger*LDi;

Psi = Lambdadagger*SigOmega(1:s,1:s)*Lambdadagger'-Di*Di';


% correct, if not positive definite.
[u,su,v]=svd(Psi);
sp = diag(max(diag(su),0));
Psi = u*sp*u';

% convert into parameter vector: 
thetai = conv_syst_par(Ai,Bi,Ci,Di,Psi);

% function to minimize: log det SigOmega(B,D,Omega,Psi) +
% tr[SigOmega(B,D,Omega,Psi)^(-1) (ZZ/T - CAi XtX CAi')]

options = optimoptions('fminunc','display','iter');
options.MaxFunctionEvaluations = 1000;

% adjust matrices for Lambdai 
LZtZL = zeros(r+n,r+n);
LZtZL(1:r,1:r) =  Lambdadagger*ZtZ(1:s,1:s)*Lambdadagger';
LZtZL((r+1):end,1:r) =  ZtZ((s+1):end,1:s)*Lambdadagger';
LZtZL(1:r,(r+1):end) =  Lambdadagger*ZtZ(1:s,(s+1):end);
LZtZL((r+1):end,(r+1):end) =  ZtZ((s+1):end,(s+1):end);

LZtX = zeros(r+n,n);
LZtX(1:r,:)= Lambdadagger*ZtX(1:s,:);
LZtX((r+1):end,:)=ZtX((s+1):end,:);

[th_opt,fval,exitflag] =  fminunc(@(x) cal_crit_EM(x,LZtZL,LZtX,XtX,r,n), thetai,options);

%th_opt = thetai;
[Ae,Be,Ce,De,Pse]= conv_par_syst(th_opt,r,n);


% renormalize.
BDtBD = [De*De',De*Be';Be*De',Be*Be'];
[U,S,V]=svd(BDtBD);
BD = [U(:,1:r)]*inv(U(1:r,1:r));
Be = BD((r+1:end),:);
De = BD(1:r,:);

Omegae = BDtBD(1:r,1:r); 

