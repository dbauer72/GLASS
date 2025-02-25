function [qlike,Fe,tres,uhat] = cal_quasi_like_aDFM_thetatau(param,Y,r,n,Pbullet)
% calculates the quasi likelihood for aDFM system for fixed Lambda matrix. 
%
% SYNTAX: [qlike,Fhat,tres,uhat] = cal_quasi_like_aDFM_thetatau(param,Y,r,n,Pbull)
%
% INPUTS:  param ... d x 1 vector of parameter values (fed into param_syst_aDFM_thetatau,
%                       see there for description)
%          Y ... NxT data matrix.
%          r ... integer; static factor dimension
%          n ... integer; system order.
%          Pbull ... indicator; if Pbull>0, state is started in stationary
%                   distribution.
%
% OUTPUT:   qlike ... real; -2/T log Gaussian likelihood.
%           Fhat  ... real matrix r x T; estimated static factors
%           tres  ... rx T; matrix of innovations for static factors.
%           uhat ...  Txq; matrix of dynamic factors. 
%
%         
% AUTHOR: dbauer, 25.2.2025


if nargin<5 % Pbullet: start with x_1 = x_bullet? 
    Pbullet = 1;
end;

qlike = 0;

% dimensions 
[N,T]=size(Y);

% extract matrices from parameters 
[th,Lambda,RN,UN,Rtilde,Utilde] = param_syst_aDFM_thetatau(param,N,r,n);

% system
Omega = Rtilde*th.Omega*Rtilde'; 
A = th.A;
B= th.B*inv(Rtilde);
C= Rtilde*th.C;
D= th.D;
Gam_zeta= eye(r)/N; 

Ddagger = inv(D'*D)*D';
K = B*Ddagger;

% estimate factors 
Fhat = (Utilde'*Y/sqrt(N))'; 
Fe = Fhat*0; 

% new switch for Pbullet == 0: PE estimation 
% --- if not minimum-phase-> penalize and project! ---
makcA = max(abs(eig(A-K*C)));
if (makcA>0.99) && (Pbullet<1)
    A = A*0.99/makcA;
    K = K*0.99/makcA;
    qlike = 10^6*(exp(makcA)-exp(0.99));
end

% add constants due to directions outside of Lambda.
% choose penalty for directions orthogonal to Lambda. 
c = 1; 
qlike = qlike + r*log(N) + c*(sum(Y(:).^2)/(T) - sum(Fhat(:).^2)/T*N);

% PE estimation for factors.
% This is a shortcut which is not unproblematic: Assumes invertibility of
% A-B*Ddagger*C and uses this for the inverse filter!
    
if Pbullet<0 
    xh = ltitr(A-K*C,K,Fhat,zeros(n,1));
    tres = Fhat-xh*C';
    uhat = tres*Ddagger';
    Fe = Fe + uhat*D'; 

    Omegah = tres'*tres/T;
    qlike = qlike+ T*log(det(Omegah)) + T*r; 
    return;
end

if Pbullet>0 % if P_bull to be included: system cannot be unstable!
    maA = max(abs(eig(A)));
    if maA >0.99
        A = A*0.99/maA;
        qlike = qlike + 10^6*(exp(maA)-exp(0.99));
    end
end

% calculate P_bull
Q = K*Omega*K';
Q = (Q+Q')/2;


pb = inv( eye(n^2) - kron(A,A))*Q(:);
P0 = reshape(pb,n,n);        
P0 = (P0+P0')/2;
% --- initialize the Kalman filter ---
x0= zeros(n,1);
tres(1,:)=Fhat(1,:)-x0'*C';
Fe(1,:)= x0'*C';

Omegat= C*P0*C'+D*Omega*D'+Gam_zeta;
Kt = (A*P0*C'+K*Omega*D')*inv(Omegat);
xf = A*x0 + Kt*tres(1,:)';
Pkg1 = A*P0*A' + Q- Kt*Omegat*Kt'; % P(2|1)

qlike = qlike +  (log(det(Omegat)) + tres(1,:)*inv(Omegat)*tres(1,:)')/T;

for t=2:T % filter equations
    Omegat = C*Pkg1*C'+D*Omega*D'+Gam_zeta;
    iOm = inv(Omegat);
    tres(t,:)= Fhat(t,:) - xf'*C';
    Fe(t,:)= xf'*C';
    Kt = (A*Pkg1*C'+K*Omega*D')*inv(Omegat);
    xf = A*xf+ Kt*tres(t,:)';
    Pkg1 = (A * Pkg1 *A') +Q - (Kt*Omegat*Kt');
    Pkg1 = (Pkg1+Pkg1')/2;
    % update likelihood
    qlike = qlike + (log(det(Omegat)) + tres(t,:)*inv(Omegat)*tres(t,:)')/T;
end

uhat = tres*Ddagger';
Fe = Fe + uhat*D'; 
