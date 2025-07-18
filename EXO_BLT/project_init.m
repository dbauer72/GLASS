function th = project_init(thi,si,sist,ni)
% project_init assumes exogeneity and projects the joint system onto 
% the GLASS regional model.
% 
% SYNTAX: th = project_init(thi,si,sist,ci,cist);
% 
% INPUTS:  thi ... initial estimate of joint system
%          si  ... integer; dimension of y_{i,t}.
%          sist ... integer; dimension of star var y_{i,t}^*
%          ni  ... integer; dimension of state of conditional model
%  
% OUTPUT: th ... theta structure 
%
% AUTHOR: dbauer, 11.7.2025.

A = thi.A;
C = thi.C;
K = thi.K;
Om = thi.Omega; 

Abar = A-K*C;
n = size(A,1);

% derive D. 
CO = chol(Om)';
D= CO((sist+1):end,1:sist); 

Df = eye(si+sist);
Df((sist+1):end,1:sist)=-D; 

% convert C to block diagonal. 
C = Df*C;

% calculate Hankel matrices 
%Hst = my_hank(Abar,K(:,1:sist),C(1:sist,:),nist,nist);
% Hankel matrix corresponding to the inverse system: 
Hi = my_hank(Abar,K,C((sist+1):end,:),5*ni,5*ni);

% get realisation of the Hankel matrix.
[U,S,V]=svd(Hi);
Of = U(:,1:ni);
Cp = Of'*Hi;

% realize system 
C = Of(1:si,:);
Abar = Of(1:(end-si),:)\Of((si+1):end,:);
B = Cp(:,1:sist);
K = Cp(:,sist+[1:si]);
A = Abar+ K*C;


%[Jit,Bit,Cit,Kit] = calculate_normed_system_from_OfCp(Of,Cp,D,si,sist,ni,nist,ci,cist);

% fill in the result. 
th = theta_urs();
th.which ='SS';
th.A = A;
th.B = B;
th.C = C;
th.D = D;
th.K = K;
th.Omega = Df((sist+1):end,:)*Om*Df((sist+1):end,:)';





