function [th,RN,UN,LambdaL] = norm_aDFM_Utilde(th,Lambda)
% norm_aDFM normalises the aDFM system given by the loading matrix Lambda
% and the (tall) state space system th. 
%
% SYNTAX:  [th,RN,Utilde,Lambda] = norm_aDFM_Utilde(th,Lambda)
%
% INPUTS:   th ... theta structure for system
%           Lambda ... Nxr loading matrix. 
%
% OUTPUTS:   th ... theta structure for system
%            RN ... rxr p.l.t part of the loading matrix. 
%            Utilde ... Nxr matrix parameterizing the column space of the
%                   loading matrix. 
%
%
% REMARK: model: y_t = Lambda F_t + chi_t; 
%  F_t = Cx_t + Du_t, x_{t+1} = A x_t + Bu_t. 
% 
% Normalisation: + Lambda'Lambda/T = I_r/N + RN*RN'.
%                + Lambda = [I;UN*RN*sqrt(N)], UN'*UN=I_r.  
%                + V(u_t) = Omega.
%                + D =I_r. 
%                + (A,B,C) in echelon form. 
% AUTHOR: dbauer, 24.2.2025

[N,r] = size(Lambda);

A = th.A;
B = th.B;
C = th.C;
D = th.D; 
Omega = th.Omega; 

[n,q] = size(B);

% normalize Lambda:
TrafoL = Lambda(1:r,1:r);
LambdaL = Lambda*inv(TrafoL);
indrev = [r:-1:1];
[Q,R] = qr(LambdaL((r+1):N,indrev));
UN = Q(:,indrev)*sqrt(N);
RN = R(indrev,indrev)/sqrt(N);
RtR = eye(r)/N + RN'*RN;
Rtilde =chol(RtR);
Utilde = LambdaL*inv(Rtilde); 

C = TrafoL*C;
Omega = TrafoL*Omega*TrafoL';

% convert state space system to echelon canonical form. 
% normalize system 
On = zeros(r*n,n);
On(1:r,:)=C; 
for j=1:(n-1)
    On(j*r+[1:r],:) = C*A^j;
end

Trafo = On(1:n,1:n);
iTrafo = inv(Trafo);

th.C = C*iTrafo;
th.A = Trafo*A*iTrafo;
th.B = Trafo*B*inv(TrafoL);
th.D = D;
th.Omega = Omega; 

