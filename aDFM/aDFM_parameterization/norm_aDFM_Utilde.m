function [th,RN,UN,LambdaL] = norm_aDFM_Utilde(th,Lambda,c)
% norm_aDFM normalises the aDFM system given by the loading matrix Lambda
% and the (tall) state space system th. 
%
% SYNTAX:  [th,RN,Utilde,Lambda] = norm_aDFM_Utilde(th,Lambda)
%
% INPUTS:   th ... theta structure for system
%           Lambda ... Nxr loading matrix. 
%           c  ... integer; number of common trends in common factors 
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
%                + (A,B,C) in Bauer Wagner canonical form. 
% AUTHOR: dbauer, 10.3.2025

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
B= B*inv(TrafoL); 

% convert state space system to echelon canonical form. 
% normalize system 
if nargin<3
    c = 0;
end;

% transform to extract common trends part. 
if (c>0)
    [V,dA] = eig(A-eye(n));
    dd = diag(dA);
    [~,I]=sort(abs(dd));
    V= V(:,I);
    A = inv(V)*A*V; 
    C = C*V;
    B = inv(V)*B; 
    % convert C_1B_1 into canonical form 
    C1B1 = real(C(:,1:c)*B(1:c,:));
    [Q,R]=qr(C1B1);
    C(:,1:c)=Q(:,1:c);
    B(1:c,:)=R(1:c,:);

    indbull= (c+1):n; 
    Abull = A(indbull,indbull);
    Bbull = B(indbull,:);
    Cbull = C(:,indbull); 
else
    Abull = A;
    Bbull = B;
    Cbull = C; 
    indbull = 1:n; 
end

On = zeros(r*(n-c),(n-c));
On(1:r,:)=Cbull; 
for j=1:(n-1)
    On(j*r+[1:r],:) = Cbull*Abull^j;
end

Trafo = On(1:(n-c),1:(n-c));
iTrafo = inv(Trafo);

tCbull = real(Cbull*iTrafo);
tAbull = real(Trafo*Abull*iTrafo);
tBbull = real(Trafo*Bbull);

% fill in new system 
C(:,indbull)= tCbull;
A(indbull,indbull)=tAbull;
B(indbull,:)=tBbull;

if (c>0)
    A(1:c,1:c) = eye(c);
    A(1:c,indbull)=0;
    A(indbull,1:c)=0; 
end

th.C = C;
th.A= A;
th.B = B; 
th.D = D;
th.Omega = Omega; 

