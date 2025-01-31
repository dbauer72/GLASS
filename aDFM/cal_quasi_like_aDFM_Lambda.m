function [qlike,Fe,tres,uhat] = cal_quasi_like_aDFM_Lambda(param,Lambda,Y,r,q,n,Pbullet)
% calculates the quasi likelihood for aDFM system for fixed Lambda matrix. 
%
% SYNTAX: [qlike,Fhat,Fres] = cal_quasi_like_aDFM_Lambda(param,Y,r,q,n,Pbull)
%
% INPUTS:  param ... d x 1 vector of parameter values (fed into param_syst_aDFM_Lambda,
%                       see there for description)
%          Y ... NxT data matrix.
%          r ... integer; static factor dimension
%          q ... integer; dynamic factor dimension
%          n ... integer; system order.
%          Pbull ... indicator; if Pbull>0, state is started in stationary
%                   distribution.
%
% OUTPUT:   qlike ... real; -2/T log Gaussian likelihood.
%           Fhat  ... real matrix r x T; estimated static factors
%           Fres  ... qx T; matrix of estimated dynamic factors. 
%
%         
% AUTHOR: dbauer, 24.7.2024


if nargin<6 % Pbullet: start with x_1 = x_bullet? 
    Pbullet = 0;
end;

qlike = 0;

% dimensions 
[N,T]=size(Y);

% extract matrices from parameters 
[th,Gam_zeta] = param_syst_aDFM_Lambda(param,N,r,q,n);
Omega = eye(q); % normalisation of variance of dynamic factors to identity matrix.
A = th.A;
B= th.B;
C= th.C;
D= th.D;

Ddagger = inv(D'*D)*D';
K = B*Ddagger;

% estimate factors 
Fhat = (Lambda'*Y/N)'; 
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
Q = B*B';
Q = (Q+Q')/2;


pb = inv( eye(n^2) - kron(A,A))*Q(:);
P0 = reshape(pb,n,n);        
P0 = (P0+P0')/2;
% --- initialize the Kalman filter ---
x0= zeros(n,1);
tres(1,:)=Fhat(1,:)-x0'*C';
Fe(1,:)= x0'*C';

Omegat= C*P0*C'+D*Omega*D'+Gam_zeta;
Kt = (A*P0*C'+B*Omega*D')*inv(Omegat);
xf = A*x0 + Kt*tres(1,:)';
Pkg1 = A*P0*A' + Q- Kt*Omegat*Kt'; % P(2|1)

qlike = qlike +  (log(det(Omegat)) + tres(1,:)*inv(Omegat)*tres(1,:)')/T;

for t=2:T % filter equations
    Omegat = C*Pkg1*C'+D*Omega*D'+Gam_zeta;
    iOm = inv(Omegat);
    tres(t,:)= Fhat(t,:) - xf'*C';
    Fe(t,:)= xf'*C';
    Kt = (A*Pkg1*C'+B*Omega*D')*inv(Omegat);
    xf = A*xf+ Kt*tres(t,:)';
    Pkg1 = (A * Pkg1 *A') +Q - (Kt*Omegat*Kt');
    Pkg1 = (Pkg1+Pkg1')/2;
    % update likelihood
    qlike = qlike + (log(det(Omegat)) + tres(t,:)*inv(Omegat)*tres(t,:)')/T;
end

uhat = tres*Ddagger';
Fe = Fe + uhat*D'; 
