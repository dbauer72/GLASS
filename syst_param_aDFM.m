function param = syst_param_aDFM(th,Lambda)
% syst_param_aDFM converts a normalized system into a parameter vector.
%
% SYNTAX:  param = syst_param_aDFM(th,Lambda)
%
% INPUTS:   th ... theta structure for system
%           Lambda ... Nxr loading matrix. 
%
% OUTPUTS: param ... parameter vector. 
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

[N,r]=size(Lambda);

% parameters for Lambda
parLam = ortho2par_plt(Lambda/sqrt(N));

% now extract parameters from system
param = [parLam(:)];
[n,q] = size(th.B);
% parameters for C
if (r>n) % only if r is larger than n the matrix C contains parameters. 
    pc = th.C((n+1):r,:);
    param = [param;pc(:)]; 
end

% parameters for A
if (n <= r)
    param = [param;th.A(:)];
else
    pa = th.A((n-r+1):end,:);
    param = [param;pa(:)];
end

% parameters for B
param = [param;th.B(:)];

% parameters for D
for j=1:q
    param = [param;th.D(j:end,j)];
end


