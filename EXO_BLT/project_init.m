function th = project_init(thi,si,sist,ni,nist,ci,cist);
% project_init assumes exogeneity and projects the joint system onto 
% the GLASS regional model.
% 
% SYNTAX: th = project_init(thi,si,sist,ci,cist);
% 
% INPUTS:  thi ... initial estimate of joint system
%          si  ... integer; dimension of y_{i,t}.
%          sist ... integer; dimension of star var y_{i,t}^*
%          ci  ... integer; number of common trends in y_{i,t}
%          cist ... integer; number of common trends in y_{i,t}^*.
%
%  
% OUTPUT: th ... theta structure 
%
% AUTHOR: dbauer, 11.7.2025.

A = thi.A;
C = thi.C;
K = thi.K;
Om = thi.Omega; 

Abar = A-K*C;
n = ni+nist;

% derive D. 
CO = chol(Om)';
D= CO((sist+1):end,1:sist); 

Df = eye(si+sist);
Df((sist+1):end,1:sist)=-D; 

% convert C to block diagonal. 
C = Df*C;

% calculate Hankel matrices 
Hst = my_hank(Abar,K(:,1:sist),C(1:sist,:),nist,nist);
Hi = my_hank(Abar,K,C((sist+1):end,:),2*ni,2*ni);

% get realisation of the Hankel matrix.
[U,S,V]=svd(Hi);
Of = U(:,1:ni);
Cp = Of'*Hi;
Ci = U(1:si,1:ni);
Ki = Cp(:,(sist+1):(sist+si)); 
Bi = Cp(:,1:sist) + Ki*D;

Abari = Of(1:(end-si),:)\Of((si+1):end,:);
Ai = Abari+ Ki*Ci;

% project to get number of common trends right. 
[v,d]= eig(Ai);
ds = abs(diag(d));
[dsort,I] = sort(ds);
[dsort,I]=sort(ds,'descend');

% calculate transformation and normalisation 
Trafo = v(:,I);
Jit = inv(Trafo)*Ai*Trafo;

% unit eigenvalue part
Jit(1:ci,1:ci)=eye(ci);
Jit(1:ci,(ci+1):end) = 0;
Jit((ci+1):end,1:ci) = 0;
Cit = Ci*Trafo;
Bit = inv(Trafo)*Bi;
Kit = inv(Trafo)*Ki;

Cit(:,1:ci)=real(Cit(:,1:ci));
Bit(1:ci,:)=real(Bit(1:ci,:));
Kit(1:ci,:)=real(Kit(1:ci,:));

[Q,R]= qr(Cit(:,1:ci));
Cit(:,1:ci)=Q(:,1:ci);
Bit(1:ci,:)=R(1:ci,1:ci)*Bit(1:ci,:);
Kit(1:ci,:)=R(1:ci,1:ci)*Kit(1:ci,:);

% stable part
ind = (ci+1):ni;

% extract matrices
Abull = Jit(ind,ind);
Cbull = Cit(:,ind);
Kbull = [Bit(ind,:),Kit(ind,:)];

% calculate impulse response 
IR =zeros(si,sist+si,2*ni);
CA = Cbull;
for j=1:(2*ni)
    IR(:,:,j)=real(CA*Kbull);
    CA = CA*Abull;
end

% calculate Hankel matrix
H = zeros(si*ni,(ni+1)*(si+sist));
for a=1:ni
    for b=1:(ni+1)
        H((a-1)*si+[1:si],(b-1)*(sist+si)+[1:(sist+si)])= squeeze(IR(:,:,a+b-1));
    end
end

[u,s,v]= svd(H);
nbull= ni- ci;
Of = u(:,1:nbull);
Cp = Of'*H;
Cbull = Of(1:si,:);
Abull = Of(1:(end-si),:)\Of((si+1):end,:);
Bbull = Cp(:,1:sist);
Kbull = Cp(:,sist+[1:si]);

% fill in the matrices 
Jit(ind,ind)= Abull;
Cit(:,ind) = Cbull;
Bit(ind,:)=Bbull;
Kit(ind,:) = Kbull;

% restriction to I(1)
Bit(1:ci,1:cist)=0;

Bit = Bit + Kit*D;

% fill in the result. 
th = theta_urs();
th.which ='SS';
th.A = Jit;
th.B = Bit;
th.C = Cit;
th.D = D;
th.K = Kit;
th.Omega = Df((sist+1):end,:)*Om*Df((sist+1):end,:)';





