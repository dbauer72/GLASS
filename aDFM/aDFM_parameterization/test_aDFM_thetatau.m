% script to test functions 
clear all

% main integers: 
T = 100;
N =100; 
q = 2;
n=8;
r=q; 

% idiosynchratic part 
ths(1) = theta_urs();
ths(1).A = 1.9*(rand(1)-.5);
ths(1).B = 1;
ths(1).C = ths(1).A;
ths(1).D = 1;
ths(1).Omega = .5; 

for j=2:N
    ths(j) = theta_urs();
    ths(j).A = 1.9*(rand(1)-.5);
    ths(j).B = 1;
    ths(j).C = ths(j).A;
    ths(j).D=1;
    ths(j).Omega = .5; 
end

% common factor part. 
Lambda = par2ortho_plt(rand(N*r),N,r)*sqrt(N);
Lambda = Lambda(N:-1:1,:);
[Q,R]=qr(Lambda(1:r,:)');
Lambda= Lambda*Q*diag(sign(diag(R)));

th = theta_urs();

th.C = randn(r,n);

th.A = diag(1.6*(rand(n,1)-.5));
th.B = par2ortho_plt(rand(N),n,q);
th.D = [rand(r-q,q);eye(q)]; 
th.Omega = eye(q);

% make sure system is invertible 
ev = abs(eig(th.A-th.B*th.C));
while (max(ev)>0.99)
    th.B = th.B/2;
    ev = abs(eig(th.A-th.B*th.C));
end

[th,RN,UN,Lambda]=  norm_aDFM_Utilde(th,Lambda);
th_chi= th;
th_chi.C = Lambda*th.C;
th_chi.D = Lambda*th.D; 

% simulate the process
[y,chi,u,x,e,F]= simu_GDFM(T,ths,th,Lambda);


[resulti,thei,thii,llei] = est_aDFM_Utilde(y,r,n);

% estimate with true Lambda
[result,the] = est_aDFM_Utilde(y,r,n,UN,1,thei);

% estimate including tau_U. 
[resultf,thef] = est_aDFM_thetatau(y,r,n,1,thei);

% check, how good the result is: 
[norm(Lambda- resulti.Lambda),norm(Lambda-result.Lambda),norm(Lambda-resultf.Lambda)]


% transfer function estimate
th.K = th.B*0;
M = 100;
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
norm_diff= zeros(10,4);
for j=1:(M-1)
    norm_diff(j,1)= norm(squeeze(IF(:,:,j+1)-IFei(:,:,j+1)));
    norm_diff(j,2)= norm(squeeze(IF(:,:,j+1)-IFi(:,:,j+1)));
    norm_diff(j,3)= norm(squeeze(IF(:,:,j+1)-IFe(:,:,j+1)));
    norm_diff(j,4)= norm(squeeze(IF(:,:,j+1)-IFef(:,:,j+1)));
end
norm_diff




