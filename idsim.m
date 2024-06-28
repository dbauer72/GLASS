function [y,x,u] = idsim(T,th);
% idsim simulates a state space process using the th theta structure and
% Gaussian white noise with variance th.Omega.
%
%  SYNTAX:    y = idsim(T,th);
%
%  INPUT:   T ... integer; sample size.
%           th ...theta structure;
%
%  OUTPUT:  y Txs real matrix of simulated data.
%
%  REMARK: Uses a system in innovation form 
%  y_t = Cx_t + Du_t, 
%  x_{t+1} = Ax_t + B u_t. 
% 
% u_t is Gaussian white noise with mean zero and variane th.Omega. 
% Dimension of u_t and y_t can be different handling the singular case. 
%
%  AUTHOR: dbauer 27.1.2023. 

% dimensions
q = size(th.B,2);
s = size(th.C,1); 
n = size(th.A,1); 

if (size(th.D,1) ~= s)
    error('Output dimensions do not match!');
end
if (size(th.D,2) ~= q)
    error('Input dimensions do not match!');
end

% inputs 
u = randn(q,T);
% take th.Omega into account 
LOm = chol(th.Omega)'; 
u = LOm*u;

x = zeros(n,T+1);
y = zeros(s,T); 

% filter 
for t=1:T
    y(:,t)= th.C*x(:,t) + th.D*u(:,t);
    x(:,t+1) = th.A*x(:,t) + th.B*u(:,t);
end

y = y';
