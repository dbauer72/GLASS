function crit = cal_crit_EM(theta,ZZ,ZX,XX,s,n);
%
% 

% lower bound on eigenvalues. 
tol = 0.000001; 
eps = 0.000000001;
% return value
crit = 0; 

% generate matrices from theta vector 
[A,B,C,D,Psi]= conv_par_syst(theta,s,n);

% check stability 
[ev] = eig(A);
dev = abs(ev); 
mdev = max(dev); 
if (mdev>.99)
    crit = crit + 10^6*(exp(mdev)-1);
end


% regularize 
SigOmega = [D*D'+Psi, D*B';B*D',B*B'+ eps*eye(n)];

% extract the dominant s evs of SigOmega:
%[U,S,V]=svd(SigOmega);
%BDt = U(:,1:s);
%
%Omega = BDt'*SigOmega*BDt;

% calculate innovation variance estimate
AC = [C;A];
ZZmAC = ZZ - AC*ZX' - ZX*AC' + AC*XX*AC'; 

% regularize, if ZZmAC has negative evs.
[eve,ev] = eig(ZZmAC);
dev = diag(ev);
if min(dev)<0
    dev = max(dev,eps);
    ZZmAC = eve*diag(dev)*eve'; 
end


% calculate criterion function
tildeI_y = [eye(s);zeros(n,s)];
tildeI_x = [zeros(s,n);eye(n)];

crit = crit + log(det(tildeI_y'*SigOmega*tildeI_y)) + trace(inv(tildeI_y'*SigOmega*tildeI_y)*tildeI_y'*ZZmAC*tildeI_y);
crit = crit + log(det(tildeI_x'*SigOmega*tildeI_x)) + trace(inv(tildeI_x'*SigOmega*tildeI_x)*tildeI_x'*ZZmAC*tildeI_x);


