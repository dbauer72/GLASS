function par = th2param_X(th,c,norm);
% converts th structure into parameter vector
% assumes echelon form for stationary system, heading cxc system to
% correspond to integrated part. Includes exogenous inputs part B,D.
%
% SYNTAX:  par = th2param_X(th,c);
%
% INPUT:   th ... theta structure
%          c  ... number of common trends
%          norm ... indicator: shall the system first be normalized to
%          canform? 
%
% OUTPUT:  param ... kd x 1 real vector containing the parameters.
%
% REMARK: currently restricted to stationary case!! c and restrict are
% inactive. 
%
% AUTHOR: dbauer, 27.4.2023.

A = th.A;
B = th.B;
C = th.C; 
D = th.D;
K = th.K;

[n,s] = size(K);
m = size(B,2);

%nbull = n-c; 

% --- check if normalization necessary? ---
if (nargin < 3) | (norm==1) 
%    [V,D] = eig(th.A);
%    [~,I]=sort(abs(eig(th.A)-1));
%    Aci = inv(V(:,I))*th.A*(V(:,I));
%    Kci = inv(V(:,I))*th.K;
%    Cci = th.C*(V(:,I));
%    
%    % now standardize integrated part
%    C1 = Cci(:,1:c);
%    K1 = Kci(1:c,:);
%    
%    [Q,R] = qr(C1*K1);
%    C1 = Q(:,1:c);
%    K1 = R(1:c,:);
%    
%    % standardize stable part
%    Abull = Aci(c+1:end,c+1:end);
%    Kbull = Kci(c+1:end,:);
%    Cbull = Cci(:,c+1:end);
%    
%    [~,Ab,Cb,Kb] = ss2ech_n(Abull',Cbull',Kbull');
%
%    % fill in 
%    A = Aci;
%    A(c+1:end,c+1:end)=Ab';
%    K = [K1;Kb'];
%    C = [C1,Cb'];
% new code for stationary ARMAX systems 
[~,A,BK,C] = ss2ech_n(A,[B,K],C);
end

% --- convert to real, if necessary ----
A= real(A);
C = real(C);
BK = real(BK);

%% --- partition in C1 and Cbull ---
%C1 = C(:,1:c);
%Cbull = C(:,c+1:end);
%
%% --- partition in K1 and Kbull ---
%K1 = K(1:c,:);
%Kbull = K(c+1:end,:);

%% --- partition into A1 and Abull ---
%Abull = A(c+1:end,c+1:end);
%
% --- convert to parameters ---
%param1 = mat2param(C1,K1);
param = mat2param_bull_X(A,BK,C);

% --- handle D ---
nD = sum(D(:).^2);
if nD>0 % D has been estimated -> include parameters. 
    par = [param(:);D(:)];
else
    par = param; 
end
par = real(par);



