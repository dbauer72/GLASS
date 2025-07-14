function [par,parI1] = th2param_BLT(th,index,norm);
% converts th structure into parameter vector for regional model (ML) in
% block lower triagnular (BLT) form the context of GLASS.
%
% The structure is characterized by the index = [c(i),c(i)^*]. 
% Here c(i) the number of common trends in the regional model i
%      c(i)^* the number of common trends in y_{i,t}^* potentialy spilling
%      over. 
%
% The diagonal model uses the I(1) structure. 
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
D = th.D;
B = th.B; 

Omega= th.Omega; 

[n,si] = size(K);
ci = index(1);
cist = index(2); 

nbull = n-ci; 

% --- check if normalization necessary? ---
if (nargin < 3) | (norm==1) 
    [V,D] = eig(th.A);
    [~,I]=sort(abs(eig(th.A)-1));
    Aci = inv(V(:,I))*th.A*(V(:,I));
    Kci = inv(V(:,I))*th.K;
    Cci = th.C*(V(:,I));
    Bci = inv(V(:,I))*th.B;

    % now standardize integrated part
    C1 = Cci(:,1:ci);
    K1 = Kci(1:ci,:);
    
    [Q,R] = qr(C1*K1);
    C1 = Q(:,1:ci);
    K1 = R(1:ci,:);
    
    Bci(1:ci,:) = C1'*Cci(:,1:ci)*Bci(1:ci,:);
    % standardize stable part
    Abull = Aci(ci+1:end,ci+1:end);
    Kbull = Kci(ci+1:end,:);
    Cbull = Cci(:,ci+1:end);
    
    [~,Ab,Cb,Kb,Tr] = ss2ech_n(Abull',Cbull',Kbull');

    % fill in 
    A = Aci;
    A(ci+1:end,ci+1:end)=Ab';
    K = [K1;Kb'];
    C = [C1,Cb'];
    B = [Bci(1:ci,:);Tr'*Bci((ci+1):end,:)];
end

% --- convert to real, if necessary ----
A= real(A);
C = real(C);
K = real(K);
B = real(B);
D= real(D);

% --- partition in C1 and Cbull ---
C1 = C(:,1:ci);
Cbull = C(:,ci+1:end);

% --- partition in K1 and Kbull ---
K1 = K(1:ci,:);
Kbull = K(ci+1:end,:);

% --- partition into A1 and Abull ---
Abull = A(ci+1:end,ci+1:end);

% --- convert to parameters ---
param1 = mat2param(C1,K1);
parambull = mat2param_bull(Abull',Cbull',Kbull');

par = [param1(:);parambull(:)];
par = real(par);

parI1 = par; 
% --- next remaining parameters in D and B ---
par = [par(:);th.D(:)];
B(1:ci,2:cist) = 0;
B1 = B(1:ci,cist+1:end);
B2 = B((ci+1):end,:);

par = [par;B1(:);B2(:)];

% --- last Omega ---
 paro = extr_lowtri(Omega);
 par = [par(:);paro(:)];


