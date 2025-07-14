function [y,chi,u,x,e,F]= simu_GDFM(T,ths,th_chi,Lambda)
%  function to simulate a common factor model. 
%
%  SYNTAX:  [y,chi,u]= simu_GDFM(T,ths,q,th_chi);
%
%  INPUTS:  T ... integer; sample size
%          ths ... vector of theta structures of length n.
%          th_chi ... theta structure for common factors.
%
%  OUTPUT:  y ... NxT real matrix; observations.
%          chi ... NxT real matrix; common factor part.
%          u ... NxT real matrix; y = chi + u idiosynchratic component. 
%
%  REMARK: It is assumed that the dimensions add up. All noise sequences
%  are drawn normally distributed independent of each other. 
%
%  AUTHOR: dbauer, 27.1.2023

% idiosynchratic component 
u = zeros(0,T); 
N = length(ths); 
for in = 1:N
    th = ths(in);
    ui = idsim(T,th); 
    u = [u;ui'];
end

N = size(u,1); % cross sectional dimension. 
% common component
C = th_chi.C;
nC = size(C,1);

if nargin<4
    Lambda = randn(N,nC); 
end

%th_chi.C = Lambda * th_chi.C;
%th_chi.D = Lambda * th_chi.D;


[F,x,e] = idsim(T,th_chi);

chi = Lambda*F';
% add the two parts
y = chi + u; 

