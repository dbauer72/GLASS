function [y,chi,u]= simu_GDFM(T,ths,th_chi);
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
n = length(ths); 
for in = 1:n
    th = ths(n);
    ui = idsim(T,th); 
    u = [u;ui'];
end

N = size(u,1); % cross sectional dimension. 
% common component
C = th_chi.C;
nC = size(C,1);

if (N>nC)
   C(end+1:N,:)=reshape(randn((N-nC)*T),N-nC,T);
end
if (N<nC)
    C = C(1:N,:);
end;

th_chi.C = C;
chi = idsim(T,th_chi)';

% add the two parts
y = chi + u; 

