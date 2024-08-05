function par = th2param_RM(th,index,norm);
% converts th structure into parameter vector for regional model (ML) in
% the context of GLASS.
%
% The structure is characterized by the index = [s(i),c(i),c(i)^*]. 
% Here s(i) denotes the dimension of the regional variables, 
%      c(i) the number of common trends in the regional model i
%      c(i)^* the rank of K_{i,c}^* indicating the common trends spilling
%      over. 
%
% Assumes echelon form for stationary system, heading cxc system to
% correspond to integrated part.
%
% SYNTAX:  par = th2param_RM(th,index,norm);
%
% INPUT:   th ... theta structure
%          index =[s(i),c(i),c(i)^*] ... vector of indices for i-th regional model.  
%          norm ... indicator: shall the system first be normalized to
%          canform? 
%
% OUTPUT:  param ... kd x 1 real vector containing the parameters.
%
% AUTHOR: dbauer, 2.8.2024. 

A = th.A;
K = th.K;
C = th.C; 
Omega= th.Omega; 

[n,s] = size(K);
si = index(1);
ci = index(2);
sist = index(3); 

nbull = n-ci; 

% --- check if normalization necessary? ---
if (nargin < 3) | (norm==1) 
    [V,D] = eig(th.A);
    [~,I]=sort(abs(eig(th.A)-1));
    Aci = inv(V(:,I))*th.A*(V(:,I));
    Kci = inv(V(:,I))*th.K;
    Cci = th.C*(V(:,I));
    
    % now standardize integrated part
    C1 = Cci(:,1:ci);
    K1 = Kci(1:ci,:);
    
    [Q,R] = qr(C1*K1);
    C1 = Q(:,1:ci);
    K1 = R(1:ci,:);

    % bring C1 to p.l.t. structure. 
    [Q,R] = qr(C1');
    C1 = R(1:ci,:)';
    K1 = Q(1:ci,1:ci)'*K1;
    
    % standardize stable part
    Abull = Aci(ci+1:end,ci+1:end);
    Kbull = Kci(ci+1:end,:);
    Cbull = Cci(:,ci+1:end);
    
    [~,Ab,Cb,Kb] = ss2ech_n(Abull',Cbull',Kbull');

    % fill in 
    A = Aci;
    A(ci+1:end,ci+1:end)=Ab';
    K = [K1;Kb'];
    C = [C1,Cb'];
end

% --- convert to real, if necessary ----
A= real(A);
C = real(C);
K = real(K);

% --- partition in C1 and Cbull ---
C1 = C(:,1:ci);
Cbull = C(:,ci+1:end);

% --- partition in K1 and Kbull ---
K1 = K(1:ci,:);
Kbull = K(ci+1:end,:);

% --- partition into A1 and Abull ---
Abull = A(ci+1:end,ci+1:end);

% --- convert to parameters ---
param1 = mat2param_RM(C1,K1,index);
parambull = mat2param_bull(Abull',Cbull',Kbull');

par = [param1(:);parambull(:)];
par = real(par);



