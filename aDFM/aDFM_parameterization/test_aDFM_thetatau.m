% script to test functions 
clear all

% main integers: 
T = 2000;
N =50; 
q = 2;
n=6;
r=4; 

% idiosynchratic part 
ths(1) = theta_urs();
ths(1).A = 1.9*(rand(1)-.5);
ths(1).B = 1;
ths(1).C = ths(1).A;
ths(1).D = 1;
ths(1).Omega = .05; 

for j=2:N
    ths(j) = theta_urs();
    ths(j).A = 1.9*(rand(1)-.5);
    ths(j).B = 1;
    ths(j).C = ths(j).A;
    ths(j).D=1;
    ths(j).Omega = .05; 
end

% common factor part. 
%Lambda = par2ortho(rand(N*r),N,r)*sqrt(N);
%Lambda = Lambda(N:-1:1,:);
Lambda = randn(N,r); 
Lambda(1,2)=0;
Lambda(2,1)=0;
Lambda(1,1)=0.5*rand(1)+.5;
Lambda(2,2)=0.5*rand(1)+.5;
[Q,R]=qr(Lambda(1:r,:)');
Lambda= Lambda*Q*diag(sign(diag(R)));

th = theta_urs();

C = randn(r,n);
A = diag(1.6*(rand(n,1)-.5));
B = par2ortho_plt(rand(N),n,q);
D = [rand(r-q,q);eye(q)]; 
th.Omega = eye(q);

% cumulate first element of F_t. 
c=1;

% make sure system is invertible 
ev = abs(eig(A-B*inv(D'*D)*D'*C));
while (max(ev)>0.99)
    B = B/2;
    ev = abs(eig(A-B*inv(D'*D)*D'*C));
end

n= n+c;
th.A = zeros(n,n); 
th.A(1:c,1:c)=eye(c);
th.A(1:c,c+1:end)=C(1:c,:);
th.A(c+1:end,c+1:end)=A; 
th.B = zeros(n,q);
th.B=[[eye(c),zeros(c,q-c)];B];
th.C = zeros(r,n);
th.C(1:c,1:c)=eye(c);
th.C(:,(c+1):end)=C; 
th.D = D; 

[th,RN,UN,Lambda]=  norm_aDFM_Utilde(th,Lambda,c);
th_chi= th;
th_chi.C = Lambda*th.C;
th_chi.D = Lambda*th.D; 

% simulate the process
[y,chi,u,x,e,F]= simu_GDFM(T,ths,th,Lambda);
 
M = 100;
IF = impulse(th,M,1);

Ts = [100,200,400,800,1600,3200];
for t=1:length(Ts)
    T = Ts(t) 
    for m=1:M
        [y,chi,u,x,e,F]= simu_GDFM(T,ths,th,Lambda);
        % subspace estimates
        [ths1,ths12,Lambdahats1] = cal_est_aDFM(y,r,n,q,0);
        [ths2,ths22,Lambdahats2] = cal_est_aDFM(y,r,n,q,1);
        [ths3] = cal_est_aDFM(y,r,n,q,-1);

        ths1.K = ths1.B*0;
        IFs1 = impulse(ths1,M,1);

        ths2.K = ths2.B*0;
        IFs2 = impulse(ths2,M,1);

        ths12.K = ths12.B*0;
        IFs12 = impulse(ths12,M,1);

        ths22.K = ths22.B*0;
        IFs22 = impulse(ths22,M,1);

        ths3.K = ths3.B*0;
        IFs3 = impulse(ths3,M,1);

        norm_diff_s((t-1)*M+m,:) = [norm(Lambda-Lambdahats1),norm(IFs1-IF,'fro'),norm(IFs2-IF,'fro'),norm(IFs12-IF,'fro'),norm(IFs22-IF,'fro'),norm(IFs3-IF,'fro')];
    end
    me_no(t,:)= mean(norm_diff_s((t-1)*M+[1:M],:));
end;

mean(norm_diff_s)
hist(norm_diff_s)
me_no
[me_no(:,2),me_no(:,6)]./(me_no(:,3)*[1,1])

%
c=1;
[resulti,thei,thii,llei] = est_aDFM_Utilde(y,r,n,c);

% estimate with true Lambda
[result,the] = est_aDFM_Utilde(y,r,n,c,UN,1,thei);

% estimate including tau_U. 
[resultf,thef] = est_aDFM_thetatau(y,r,n,c,1,thei);

% check, how good the result is: 
[norm(Lambda- resulti.Lambda),norm(Lambda-result.Lambda),norm(Lambda-resultf.Lambda)]/norm(Lambda)


% transfer function estimate
th.K = th.B*0;
M = 20;
IF = impulse(th,M);

thei.K = thei.B*0;
IFei = impulse(thei,M);

the.K = the.B*0;
IFe = impulse(the,M);

thii.K = thii.B*0;
IFi = impulse(thii,M);

thef.K = thef.B*0;
IFef = impulse(thef,M);


[norm(IF-IFei,'fro'),norm(IF-IFi,'fro'),norm(IF-IFe,'fro'),norm(IF-IFef,'fro')]
norm_diff= zeros(M-1,4);
for j=1:(M-1)
    norm_diff(j,1)= norm(squeeze(IF(:,:,j+1)-IFei(:,:,j+1)))/norm(squeeze(IF(:,:,j+1)));
    norm_diff(j,2)= norm(squeeze(IF(:,:,j+1)-IFi(:,:,j+1)))/norm(squeeze(IF(:,:,j+1)));
    norm_diff(j,3)= norm(squeeze(IF(:,:,j+1)-IFe(:,:,j+1)))/norm(squeeze(IF(:,:,j+1)));
    norm_diff(j,4)= norm(squeeze(IF(:,:,j+1)-IFef(:,:,j+1)))/norm(squeeze(IF(:,:,j+1)));
end
norm_diff

% same with c=0;
c=0;
[resulti0,thei0,thii0,llei0] = est_aDFM_Utilde(y,r,n,c);

% estimate with true Lambda
[result0,the0] = est_aDFM_Utilde(y,r,n,c,UN,-1,thei);

% estimate including tau_U. 
[resultf0,thef0] = est_aDFM_thetatau(y,r,n,c,-1,thei);

% check, how good the result is: 
[norm(Lambda- resulti0.Lambda),norm(Lambda-result0.Lambda),norm(Lambda-resultf0.Lambda)]/norm(Lambda)


% transfer function estimate
th.K = th.B*0;
%M = 20;
IF = impulse(th,M);

thei0.K = thei0.B*0;
IFei0 = impulse(thei0,M);

the0.K = the0.B*0;
IFe0 = impulse(the0,M);

thii0.K = thii0.B*0;
IFi0 = impulse(thii0,M);

thef0.K = thef0.B*0;
IFef0 = impulse(thef0,M);


[norm(IF-IFei0,'fro'),norm(IF-IFi0,'fro'),norm(IF-IFe0,'fro'),norm(IF-IFef0,'fro')]
norm_diff0= zeros(M-1,4);
for j=1:(M-1)
    norm_diff0(j,1)= norm(squeeze(IF(:,:,j+1)-IFei0(:,:,j+1)))/norm(squeeze(IF(:,:,j+1)));
    norm_diff0(j,2)= norm(squeeze(IF(:,:,j+1)-IFi0(:,:,j+1)))/norm(squeeze(IF(:,:,j+1)));
    norm_diff0(j,3)= norm(squeeze(IF(:,:,j+1)-IFe0(:,:,j+1)))/norm(squeeze(IF(:,:,j+1)));
    norm_diff0(j,4)= norm(squeeze(IF(:,:,j+1)-IFef0(:,:,j+1)))/norm(squeeze(IF(:,:,j+1)));
end
[norm_diff,norm_diff0]



