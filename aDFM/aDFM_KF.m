function [xt,Pt,Kt,Omegat,et,xs,Ps,xtt,Pst]= aDFM_KF(y,Lambda,A,B,C,D,Omega,Psi);
% aDFM_KF runs the Kalman filter for the aDFM system
% 
% SYNTAX: [xt,Pt,Kt,et,xs]= aDFM_KF(y,Lambda,A,B,C,D,Omega,Psi);
% 
% INPUT: y ... NxT matrix of observations 
%        Lambda ... Nxr loading matrix; Lambda'Lambda = I_r is assumed.
%        (A,B,C,D) ... state space system
%        Omega ... r x r variance of dynamic factors
%        Psi ... rxr variance of idiosyncratic process in the direction of
%        Lambda.
%
% OUTPUT: xt ... n x T matrix of filtered states
%         Pt ... n x n x T array of state estimation variance
%         Kt ... n x r x T matrix of Kalman gains
%         Omegat ... rxrxT matrix of innovation matrices. 
%         et ... N x T matrix of residuals.
%         xs ... n x T matrix of smoothed states
%
% AUTHOR: dbauer, 25.11.2024


[N,T]= size(y);
[n,r] = size(B);

if nargin<8
    Psi = eye(r)/N; 
end

if nargin<7
    Omega = eye(r);
end

% initialize objects 
xt = zeros(n,T+1);
Pt = zeros(n,n,T+1);
Kt = zeros(n,r,T);
et = zeros(r,T);
xs = xt;
xtt = xt; 

% variance matrices 
Q = B*Omega*B';
R = D*Omega*D' + Psi/N;
S = B*Omega*D'; 

% convert observations to hat F_t of smaller size
Ft = inv(Lambda'*Lambda)*Lambda'*y;
% calculate initial estimate for Pt from stationary distribution. 
Pc = reshape(inv(eye(n^2)-kron(A,A))*Q(:),n,n);

Pt(:,:,1)=Pc;
et(:,1)=Ft(:,1);
%K = (A*Pc*C' + S)*inv(C*Pc*C'+R);
% filter: forward pass 
for t=1:T   
    e(:,t) = Ft(:,t)-C*xt(:,t);
    Omegac = C*Pc*C'+R; 
    Omegat(:,:,t)= Omegac;
    K = (A*Pc*C' + S)*inv(Omegac);
    Kt(:,:,t)=K;
    xt(:,t+1)= A*xt(:,t) + K*e(:,t);
    xtt(:,t)  = xt(:,t) + Pc*C'*inv(Omegac)*e(:,t);
    Pcn = A*Pc*A' + Q - (K*Omegac*K'); 
    Pc = 0.5*(Pcn+Pcn');
    Pt(:,:,t+1) = Pc; 
end

% smoother: backward pass.
Ps = zeros(n,n); 
Psc = Ps; 
Pst = Ps;
xs(:,T+1)=xt(:,T+1);
for t= T:-1:1
    Pc = squeeze(Pt(:,:,t));
    Ptt = Pc - Pc*C'*inv(squeeze(Omegat(:,:,t)))*C*Pc;
    barA = A - squeeze(Kt(:,:,t))*C;
    Ptp1 = squeeze(Pt(:,:,t+1));
    if min(eig(Ptp1))<0.000000001
        PiPt = squeeze(Pt(:,:,t))*inv(Ptp1+0.000000001*eye(n));
    else
        PiPt = squeeze(Pt(:,:,t))*inv(Ptp1);
    end

    Ac = barA*PiPt;
    xs(:,t)= xtt(:,t) + Ac*(xs(:,t+1)-xt(:,t+1));
    Psc=  Ptt + Ac*(Psc-squeeze(Pt(:,:,t+1)))*Ac'; 
    Ps(:,:,t)= 0.5*(Psc+Psc');
    Pst(:,:,t)= Ac*squeeze(Ps(:,:,t));
end


     

