% script to test functions 
% main integers: 
T = 1000;
N = 100; 
q = 2;
n=2;
r=n+q; 

% idiosynchartic part 
ths(1) = theta_urs();
ths(1).A = 1.9*(rand(1)-.5);
ths(1).B = 1;
ths(1).C = ths(1).A;
ths(1).D = 1;
ths(1).Omega = .1; 

for j=2:N
    ths(j) = theta_urs();
    ths(j).A = 1.9*(rand(1)-.5);
    ths(j).B = 1;
    ths(j).C = ths(j).A;
    ths(j).D=1;
    ths(j).Omega = .1; 
end

% common factor part. 
Lambda = par2ortho_plt(rand(N*r),N,r)*sqrt(N);
Lambda = Lambda(N:-1:1,:);
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
[y,chi,u]= simu_GDFM(T,ths,th,Lambda);

Gam = chi*chi'/T;
[~,s,~]= svd(Gam)
diag(s)'


[resulti,thei,thii,llei] = est_aDFM_Lambda(y,r,q,n,1);
result= resulti;
for j=1:5
    Fe = result.Fhat;
    Lame = (Fe\y')';
    [the,Lambdae] = norm_aDFM(result.th,Lame);
    [result,the,thi,lle] = est_aDFM_Lambda(y,r,q,n,1,the,Lambdae);
end

% compare to full model
[resultf,thef,Lambdaf,thif,llef] = est_aDFM(y,r,q,n,1);

[resultf2,thef2,Lambdaf2,thif2,llef2] = est_aDFM(y,r,q,n,1,thei,resulti.Lambda);
[resultf3,thef3,Lambdaf3,thif3,llef3] = est_aDFM(y,r,q,n,1,thei,result.Lambda);

[llei,lle,llef3,llef2]

[norm(Lambda-resulti.Lambda),norm(Lambda-result.Lambda),norm(Lambda-resultf3.Lambda),norm(Lambda-resultf2.Lambda)]

[norm(Lambda*th.C*th.B-resulti.Lambda*thi.C*thi.B),norm(Lambda*th.C*th.B-result.Lambda*result.th.C*result.th.B),norm(Lambda*th.C*th.B-resultf3.Lambda*resultf3.th.C*resultf3.th.B),norm(Lambda*th.C*th.B-resultf2.Lambda*resultf2.th.C*resultf2.th.B)]
[norm(Lambda*th.D-resulti.Lambda*thi.D),norm(Lambda*th.D-result.Lambda*result.th.D),norm(Lambda*th.D-resultf3.Lambda*resultf3.th.D),norm(Lambda*th.D-resultf2.Lambda*resultf2.th.D)]

[th.A,thei.A,the.A,thef3.A,thef2.A]
[th.B,thei.B,the.B,thef3.B,thef2.B]
[th.C,thei.C,the.C,thef3.C,thef2.C]
[th.D,thei.D,the.D,thef3.D,thef2.D]

[norm(th.A-thef3.A),norm(th.B-thef3.B),norm(th.C-thef3.C),norm(th.D-thef3.D)]

figure
plot(resultf3.Lambda*thef3.C*thef3.B,Lambda*th.C*th.B,'.')


figure
plot(resultf3.Lambda*thef3.C*thef3.A*thef3.B,Lambda*th.C*th.A*th.B,'.')

figure
plot(resultf3.Lambda*thef3.D,Lambda*th.D,'r.')


%%% system from Barigozzi, Hallin, Luciani, Zafaroni. 
% q=1: x_{it} = a_{i}/(1-z alpha_i)u_t + zeta_{it}
T = 120;
N= 120; 

u = randn(T,1); 
zeta = randn(N,T);

for n=1:N
    a(n)= randn(1)+1;
    alphai(n)= 0.1+ 0.7*rand(1);
    chi(n,:)= filter(1,[1,-alphai(n)],u);
    sig_chi = var(chi(n,:));
    sig_zeta = var(zeta(n,:));
    zeta(n,:)=zeta(n,:)*sqrt((sig_chi/sig_zeta)*0.5);
end

y = zeta + chi;

Gam = chi*chi'/T;
[~,s,~]= svd(Gam)
diag(s)'

% estimation 
Gamy = y*y'/T; 
[~,s,~]= svd(Gam)
diag(s)'

r=6;
q= 1;
n = 5;
[resulti,thei,thii,llei] = est_aDFM_Lambda(y,r,q,n,1);
result= resulti;
for j=1:5
    Fe = result.Fhat;
    Lame = (Fe\y')';
    [the,Lambdae] = norm_aDFM(result.th,Lame);
    [result,the,thi,lle] = est_aDFM_Lambda(y,r,q,n,1,the,Lambdae);
end

% compare to full model
[resultf,thef,Lambdaf,thif,llef] = est_aDFM(y,r,q,n,1);

[resultf2,thef2,Lambdaf2,thif2,llef2] = est_aDFM(y,r,q,n,1,thei,resulti.Lambda);
[resultf3,thef3,Lambdaf3,thif3,llef3] = est_aDFM(y,r,q,n,1,thei,result.Lambda);

[llei,lle,llef3,llef2]

chiesti = resulti.Lambda*resulti.Fhat';
chieste = result.Lambda*result.Fhat';
chiestf2 = resultf2.Lambda*resultf2.Fhat';
chiestf3 = resultf3.Lambda*resultf3.Fhat';

[norm(chi-chiesti),norm(chi-chieste),norm(chi-chiestf2),norm(chi-chiestf3)]

% compare it to metric in BHLZ 
SMSEi = sum((chi(:)-chiesti(:)).^2) / sum(chi(:).^2);
SMSEe = sum((chi(:)-chieste(:)).^2) / sum(chi(:).^2);
SMSEf2 = sum((chi(:)-chiestf2(:)).^2) / sum(chi(:).^2);
SMSEf3 = sum((chi(:)-chiestf3(:)).^2) / sum(chi(:).^2);

[SMSEi,SMSEe,SMSEf2,SMSEf3]

% different values of r,n. Replicate M times  
clear all
M = 10; 
T = 120;
N= 120;
r=3;
q= 1;
nest = r-2;
h = waitbar(0,'Please wait...');
for m=1:M
    waitbar(m/M,h, sprintf('Run: %d / %d',m,M));
    % generate data
    u1 = randn(T,1);
    u2 = randn(T,1);
    
    zeta = randn(N,T);

    for n=1:N
        a(n,1)= randn(1)+1;
        alphai(n,1)= 0.1+ 0.7*rand(1);
        chi(n,:)= a(n,1)*filter(1,[1,-alphai(n,1)],u1)';
        %a(n,2)= randn(1)+1;
        %alphai(n,2)= 0.1+ 0.7*rand(1);
        %chi(n,:)= chi(n,:) + a(n,2)*filter(1,[1,-alphai(n,2)],u2)';
        sig_chi = var(chi(n,:));
        sig_zeta = var(zeta(n,:));
        zeta(n,:)=zeta(n,:)*sqrt((sig_chi/sig_zeta)*0.5);
    end

    y = zeta + chi;

    [resulti,thei,thii,llei] = est_aDFM_Lambda(y,r,q,nest,1);
    chiesti = resulti.Lambda*resulti.Fhat';
    SMSEi(m) = sum((chi(:)-chiesti(:)).^2) / sum(chi(:).^2);

    [resultf,thef,Lambdaf,thif,llef] = est_aDFM(y,r,q,nest,1);
    chiestf = resultf.Lambda*resultf.Fhat';
    SMSEf(m) = sum((chi(:)-chiestf(:)).^2) / sum(chi(:).^2);

end
close(h);

[mean(SMSEi),mean(SMSEf)]

figure
plot(resultf3.Lambda*thef3.C*thef3.B,Lambda*th.C*th.B,'.')


figure
plot(resultf3.Lambda*thef3.C*thef3.A*thef3.B,Lambda*th.C*th.A*th.B,'.')

figure
plot(resultf3.Lambda*thef3.D,Lambda*th.D,'r.')
