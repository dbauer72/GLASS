function [th,Lambda] = norm_aDFM(th,Lambda)
% norm_aDFM normalises the aDFM system given by the loading matrix Lambda
% and the (tall) state space system th. 
%
% SYNTAX:  [th,Lambda] = norm_aDFM(th,Lambda)
%
% INPUTS:   th ... theta structure for system
%           Lambda ... Nxr loading matrix. 
%
% OUTPUTS:   th ... theta structure for system
%           Lambda ... Nxr loading matrix. 
%
%
% REMARK: model: y_t = Lambda F_t + chi_t; 
%  F_t = Cx_t + Du_t, x_{t+1} = A x_t + Bu_t. 
% 
% Normalisation: + Lambda'Lambda/T = I_r.
%                + Lambda p.l.t. (heading square matrix lower triangular
%                  with positive elements on diagonal. 
%                + V(u_t)=I_q. 
%                + D is p.l.t. 
%                + (A,B,C) in echelon form. 
% AUTHOR: dbauer, 24.7.2024

[N,r] = size(Lambda);

A = th.A;
B = th.B;
C = th.C;
D = th.D; 

[n,q] = size(B);

% normalize Lambda:
[u,s,v] = svd(Lambda);

Lambda = u(:,1:r)*sqrt(N);

% p.l.t. 
[Q,R]=qr(Lambda(1:r,1:r)');
sR = diag(sign(diag(R)));
Lambda = Lambda*Q*sR; 
TrafoL = sR*Q'*s(1:r,1:r)*v(:,1:r)'/sqrt(N);


C = TrafoL*C;
D = TrafoL*D;

% normalize D
[Q,R]= qr(D(1:q,1:q)');
sR = diag(sign(diag(R)));
D = D*Q*sR;

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
th.B = Trafo*B*Q*sR;
th.D = D;

