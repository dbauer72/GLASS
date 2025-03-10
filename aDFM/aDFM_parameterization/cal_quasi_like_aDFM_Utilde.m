function [qlike,Fe,tres,uhat] = cal_quasi_like_aDFM_Utilde(param,UN,Y,r,n,c,Pbullet)
% calculates the quasi likelihood for aDFM system for fixed Lambda matrix. 
%
% SYNTAX: [qlike,Fhat,Fres] = cal_quasi_like_aDFM_Utilde(param,Utilde,Y,r,n,c,Pbull)
%
% INPUTS:  param ... d x 1 vector of parameter values (fed into param_syst_aDFM_Utilde,
%                       see there for description)
%          Utilde ... Nxr orthonormal column parameterizing column space of
%                      loading matrix.
%          Y ... NxT data matrix.
%          r ... integer; static factor dimension
%          n ... integer; system order.
%          c ... integer; number of common trends
%          Pbull ... indicator; if Pbull>0, state is started in stationary
%                   distribution.
%
% OUTPUT:   qlike ... real; -2/T log Gaussian likelihood.
%           Fhat  ... real matrix r x T; estimated static factors
%           Fres  ... qx T; matrix of estimated dynamic factors. 
%
%         
% AUTHOR: dbauer, 10.3.2025

if nargin<7 % Pbullet: start with x_1 = x_bullet? 
    Pbullet = 0;
end;

qlike = 0;

% dimensions 
[N,T]=size(Y);

% extract matrices from parameters 
[th,RN] = param_syst_aDFM_Utilde(param,r,n,c);
Lambda = [eye(r);UN*RN];
[Q,R]=qr(Lambda);
Utilde = Q(:,1:r); 
Rtilde = Utilde'*Lambda/sqrt(N);

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

% take into account possible integration
if c>0 % there is an integration going on!
    % find tilde k. via C1:
    C1 = C(:,1:c);
    Pi = C1*inv(C1'*C1)*C1';
    Cbull = [C(:,c+1:end)];
    Kbull= K(c+1:end,:);
    K1 = K(1:c,:);
    
    % calculate tilde y_t = y_t - Pi y_{t-1}
    ty = Fhat;
    ty(2:end,:)=ty(2:end,:)-Fhat(1:end-1,:)*Pi;
    
    % transformed matrices
    tilA = [zeros(c,c),-C1'*Cbull;A(c+1:end,:)];
    tilK = [K1-C1';Kbull];
    tilC = C;
else % no integration
    tilA = A;
    tilK= K;
    tilC=C;
    ty = Fhat;
end;

if Pbullet>0 % if P_bull to be included: system cannot be unstable!
    maA = max(abs(eig(tilA)));
    if maA >0.999
        tilA = tilA*0.999/maA;
        qlike = qlike + 10^6*(exp(maA)-exp(0.999));
    end
end

% new switch for Pbullet == 0: PE estimation 
% --- if not minimum-phase-> penalize and project! ---
makcA = max(abs(eig(tilA-tilK*tilC)));
if (makcA>0.99) && (Pbullet<1)
    tilA = tilA*0.99/makcA;
    tilK = tilK*0.99/makcA;
    qlike = 10^6*(exp(makcA)-exp(0.99));
end


% add constants due to directions outside of Lambda.
% choose penalty for directions orthogonal to Lambda. 
con = 1; 
qlike = qlike + r*log(N) + con*(sum(Y(:).^2)/(T) - sum(Fhat(:).^2)/T*N);

% PE estimation for factors.
% This is a shortcut which is not unproblematic: Assumes invertibility of
% A-B*Ddagger*C and uses this for the inverse filter!
    
if Pbullet<0 
    xh = ltitr(tilA-tilK*tilC,tilK,Fhat,zeros(n,1));
    tres = Fhat-xh*tilC';
    uhat = tres*Ddagger';
    Fe = Fe + uhat*D'; 

    Omegah = tres'*tres/T;
    qlike = qlike+ T*log(det(Omegah)) + T*r; 
    return;
end

% calculate P_bull
Q = tilK*Omega*tilK';
Q = (Q+Q')/2;

Qbull = Q(c+1:end,c+1:end);

pb = inv( eye((n-c)^2) - kron(tilA(c+1:end,c+1:end),tilA(c+1:end,c+1:end)))*Qbull(:);
Pbull = reshape(pb,n-c,n-c);
P0 = [0*eye(c),zeros(c,n-c);zeros(n-c,c),Pbull];


P0 = (P0+P0')/2; % P(1|0).
% --- initialize the Kalman filter ---
x0= zeros(n,1); % x(1|0)

% --- run the Kalman filter ---
tres = ty*0;
tres(1,:)=ty(1,:); % e(1)
Fe(1,:)= x0'*C';
Omegat= tilC*P0*tilC'+D*Omega*D'+Gam_zeta; % Omega(1|0)
Kt = (tilA*P0*tilC'+tilK*Omegat*D')*inv(Omegat);
xf = tilA*x0 + Kt*tres(1,:)';
Pkg1 = tilA*P0*tilA' + Q- Kt*Omegat*Kt'; % P(2|1)

ev = min(eig(Pkg1));
if ev<0
    Pkg1 = Pkg1 -eye(n)*ev; 
end

qlike = qlike +  log(det(Omegat)) + tres(1,:)*inv(Omegat)*tres(1,:)';

for t=2:T % filter equations
    Omegat = tilC*Pkg1*tilC'+D*Omega*D'+Gam_zeta;
    iOm = inv(Omegat);
    tres(t,:)= ty(t,:) - xf'*tilC';
    Fe(t,:)= xf'*tilC';
    Kt = (tilA*Pkg1*tilC'+tilK*Omega*D')*iOm;
    xf = tilA*xf+ Kt*tres(t,:)';
    Pkg1 = (tilA * Pkg1 *tilA') +Q - (Kt*Omegat*Kt');
    Pkg1 = (Pkg1+Pkg1')/2;
    % update likelihood
    qlike = qlike + log(det(Omegat)) + tres(t,:)*inv(Omegat)*tres(t,:)';
end

uhat = tres*Ddagger';
Fe = Fe + uhat*D'; 
