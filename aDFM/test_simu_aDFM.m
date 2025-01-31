% script to test functions 
% main integers: 
T = 1000;
N =20; 
q = 2;
n=2;
r=q; 

% idiosynchartic part 
ths(1) = theta_urs();
ths(1).A = 1.9*(rand(1)-.5);
ths(1).B = 1;
ths(1).C = ths(1).A;
ths(1).D = 1;
ths(1).Omega = .001; 

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
[Q,R]=qr(Lambda(1:r,:)');
Lambda= Lambda*Q*diag(sign(diag(R)));

th = theta_urs();
th.C = [eye(n);rand(r-n,n)];
th.A = diag(1.6*(rand(n,1)-.5));
th.B = par2ortho_plt(rand(N),n,q);
th.D = [rand(r-q,q);eye(q)]; 
th.Omega = eye(q);

[th,Lambda]=  norm_aDFM(th,Lambda);
th_chi= th;
th_chi.C = Lambda*th.C;
th_chi.D = Lambda*th.D; 

% simulate the process
[y,chi,u,x,e,F]= simu_GDFM(T,ths,th,Lambda);
