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
[y,chi,u,x,e,F]= simu_GDFM(T,ths,th,Lambda);

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

% alternative procedures:
% CVA for PCAs
[th_cva,rest_cva,Lambda_cva,pest,nest,nest2,s,qest,Omega] = cal_est_GDFM(y,r,n,q);
[th_cva,Lambda_cva] = norm_aDFM(th_cva,Lambda_cva);

[th_VAR,rest,Lambda_VAR,pest] = cal_est_GDFM_VAR(y,r,n,q);
[th_VAR,Lambda_VAR] = norm_aDFM(th_VAR,Lambda_VAR);

% compare to full model
[resultf,thef,Lambdaf,thif,llef] = est_aDFM(y,r,q,n,1);

[resultf2,thef2,Lambdaf2,thif2,llef2] = est_aDFM(y,r,q,n,1,thei,resulti.Lambda);
[resultf3,thef3,Lambdaf3,thif3,llef3] = est_aDFM(y,r,q,n,1,thei,result.Lambda);

[results,the,Lambdae,thi,qlike,Gam_zeta] =  aDFM_EM(y,n,r,100,-1);
[results2,the2,thi,qlike] = est_aDFM_Lambda(y,r,q,n,1,the,Lambdae,Gam_zeta)



[llei,lle,llef3,llef2,qlike]

[norm(Lambda-results.Lambda),norm(Lambda-result.Lambda),norm(Lambda-resultf3.Lambda),norm(Lambda-resultf2.Lambda)]

[norm(Lambda*th.C*th.B-results.Lambda*the.C*the.B),norm(Lambda*th.C*th.B-result.Lambda*result.th.C*result.th.B),norm(Lambda*th.C*th.B-resultf3.Lambda*resultf3.th.C*resultf3.th.B),norm(Lambda*th.C*th.B-resultf2.Lambda*resultf2.th.C*resultf2.th.B)]
[norm(Lambda*th.D-results.Lambda*the.D),norm(Lambda*th.D-result.Lambda*result.th.D),norm(Lambda*th.D-resultf3.Lambda*resultf3.th.D),norm(Lambda*th.D-resultf2.Lambda*resultf2.th.D)]

[th.A,thei.A,the.A,thef3.A,thef2.A]
[th.B,thei.B,the.B,thef3.B,thef2.B]
[th.C,thei.C,the.C,thef3.C,thef2.C]
[th.D,thei.D,the.D,thef3.D,thef2.D]
[th.Omega,thei.Omega,the.Omega,thef3.Omega,thef2.Omega]

[norm(th.A-thef3.A),norm(th.B-thef3.B),norm(th.C-thef3.C),norm(th.D-thef3.D)]

figure
plot(resultf.Lambda*thef.C*thef.B,Lambda*th.C*th.B,'.')


figure
plot(resultf3.Lambda*thef3.C*thef3.A*thef3.B,Lambda*th.C*th.A*th.B,'.')

figure
plot(resultf3.Lambda*thef3.D,Lambda*th.D,'r.')


%%% system from Barigozzi, Hallin, Luciani, Zafaroni. 
% q=1: x_{it} = a_{i}/(1-z alpha_i)u_t + zeta_{it}
clear all;
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

[results,the,Lambdae,thi,qlike,Gam_zeta] =  aDFM_EM(y,n,r,100,-1);
[results2,the2,thi,qlike] = est_aDFM_Lambda(y,r,q,n,1,the,Lambdae,Gam_zeta);


[llei,lle,llef3,llef2,qlike]

chiesti = resulti.Lambda*resulti.Fhat';
chieste = result.Lambda*result.Fhat';
chiestf2 = resultf2.Lambda*resultf2.Fhat';
chiestf3 = resultf3.Lambda*resultf3.Fhat';
chiests = results.Lambda*results.Fhat;

[norm(chi-chiesti),norm(chi-chieste),norm(chi-chiestf2),norm(chi-chiestf3), norm(chi-chiests)]

% compare it to metric in BHLZ 
SMSEi = sum((chi(:)-chiesti(:)).^2) / sum(chi(:).^2);
SMSEe = sum((chi(:)-chieste(:)).^2) / sum(chi(:).^2);
SMSEf2 = sum((chi(:)-chiestf2(:)).^2) / sum(chi(:).^2);
SMSEf3 = sum((chi(:)-chiestf3(:)).^2) / sum(chi(:).^2);
SMSEs = sum((chi(:)-chiests(:)).^2) / sum(chi(:).^2);

[SMSEi,SMSEe,SMSEf2,SMSEf3,SMSEs]

% different values of r,n. Replicate M times  
clear all
M = 10; 
T = 480;
N=480;
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
    
    %nest=2;
    %r=3
    %r=3;
    %q= 1;
    %nest = r-1;
    
    %[resulti,thei,thii,llei] = est_aDFM_Lambda(y,r,q,nest,1);
    %chiesti = resulti.Lambda*resulti.Fhat';
    %SMSEi(m,1) = sum((chi(:)-chiesti(:)).^2) / sum(chi(:).^2);

    % r= 5
    r=5;
    q= 1;
    nest = 4;
    %nest = r-1;
    
    [resulti,thei,thii,llei] = est_aDFM_Lambda(y,r,q,nest,1);
    chiesti = resulti.Lambda*resulti.Fhat';
    uhat = resulti.uhat; 
    SMSEi(m,2) = sum((chi(:)-chiesti(:)).^2) / sum(chi(:).^2);
    SMSEi(m,1) = (u1'*uhat)^2/(uhat'*uhat* u1'*u1);

    % r= 7
    %r=7;
    %q= 1;
    %nest = r-1;
    
    %[resulti,thei,thii,llei] = est_aDFM_Lambda(y,r,q,nest,1);
    %chiesti = resulti.Lambda*resulti.Fhat';
    %SMSEi(m,3) = sum((chi(:)-chiesti(:)).^2) / sum(chi(:).^2);

    %[resultf,thef,Lambdaf,thif,llef] = est_aDFM(y,r,q,nest,1);
    %chiestf = resultf.Lambda*resultf.Fhat';
    %SMSEf(m) = sum((chi(:)-chiestf(:)).^2) / sum(chi(:).^2);
end
close(h);

mean(SMSEi)

% same for q=2 
% different values of r,n. Replicate M times  
clear all
M = 1000; 
T = 480;
N=T;
h = waitbar(0,'Please wait...');
for m=1:M
    waitbar(m/M,h, sprintf('Run: %d / %d',m,M));
    % generate data
    u1 = randn(T,1);
    u2 = randn(T,1);
    u = [u1,u2]; 
    zeta = randn(N,T);

    for n=1:N
        if (n<N/2)
            alphai(n,1)=1;
            alphai(n,2)=1;
        else
            alphai(n,1)= 0.1+ 0.7*rand(1); 
            alphai(n,2)= 0.1+ 0.7*rand(1);
        end
        a(n,1)= randn(1)+1;
        
        chi(n,:)= a(n,1)*filter(1,[1,-alphai(n,1)],u1)';
        a(n,2)= randn(1)+1;
        
        chi(n,:)= chi(n,:) + a(n,2)*filter(1,[1,-alphai(n,2)],u2)';
        sig_chi = var(chi(n,:));
        sig_zeta = var(zeta(n,:));
        if (n>N/2)
            zeta(n,:)=zeta(n,:)*sqrt((sig_chi/sig_zeta)*0.5);
        end
    end

    y = zeta + chi;

    r=6;
    q= 2;
    nest = 4;
    %nest = r-1;
    
    %[resulti,thei,thii,llei] = est_aDFM_Lambda(y,r,q,nest,1);
    %chiesti = resulti.Lambda*resulti.Fhat';
    %uhat = resulti.uhat; 
    %SMSEi(m,2) = sum((chi(:)-chiesti(:)).^2) / sum(chi(:).^2);
    %SMSEi(m,1) = trace(u'*uhat*inv(uhat'*uhat)*uhat'*u)/trace(u'*u);

    %result= resulti;
    %for j=1:3
    %    Fe = result.Fhat;
    %    Lame = (Fe\y')';
    %    [the,Lambdae] = norm_aDFM(result.th,Lame);
    %    [result,the,thi,lle] = est_aDFM_Lambda(y,r,q,nest,1,the,Lambdae);
    %end

    [th,rest,Lambdahat,pest,nest,nest2,s,qest,Omega] = cal_est_GDFM(y,r,nest,q);
    Aest(m,:)= (th.A(:))';

    %chiest = result.Lambda*result.Fhat';
    %uhat = result.uhat; 
    %SMSEi(m,4) =sum((chi(:)-chiest(:)).^2) / sum(chi(:).^2);
    %SMSEi(m,3) = trace(u'*uhat*inv(uhat'*uhat)*uhat'*u)/trace(u'*u);

    %[resultf2,thef2,Lambdaf2,thif2,llef2] = est_aDFM(y,r,q,nest,1,thei,resulti.Lambda);

    %chiesti = resultf2.Lambda*resultf2.Fhat';
    %uhat = resultf2.uhat; 
    %SMSEi(m,4) = sum((chi(:)-chiesti(:)).^2) / sum(chi(:).^2);
    %SMSEi(m,3) = trace(u'*uhat*inv(uhat'*uhat)*uhat'*u)/trace(u'*u);


    % r= 7
    %r=7;
    %q= 1;
    %nest = r-1;
    
    %[resulti,thei,thii,llei] = est_aDFM_Lambda(y,r,q,nest,1);
    %chiesti = resulti.Lambda*resulti.Fhat';
    %SMSEi(m,3) = sum((chi(:)-chiesti(:)).^2) / sum(chi(:).^2);

    %[resultf,thef,Lambdaf,thif,llef] = est_aDFM(y,r,q,nest,1);
    %chiestf = resultf.Lambda*resultf.Fhat';
    %SMSEf(m) = sum((chi(:)-chiestf(:)).^2) / sum(chi(:).^2);
end
close(h);

mean(SMSEi)


figure;hold on; freq= 0:.01:2*pi;om = exp(freq*sqrt(-1));x = [real(om)',imag(om)'];plot(x(:,1),x(:,2));

for j=1:M,A= reshape(Aest(j,:),4,4);ev= eig(A);plot(real(ev),imag(ev),'x');end;

% order eigenvalues 
for j=1:M,A= reshape(Aest(j,:),4,4);ev= eig(A);evs(j,:)=ev(:);end;
for j=1:M,
    ev = evs(j,:); 
    mev = abs(ev-1);
    [mevs,I]= sort(mev,'ascend');
    evs(j,:)=evs(j,I);
    dists(j,:)= mevs;
end

dis = 0.5*(dists(:,1)+dists(:,2));
[xg,yg]=kdfft1(dis,'knorm',1024,0.001);

% compare to simulated distribution 
B = 10000;
for b=1:B
    x = cumsum(randn(500,2));
    rhomat = x(1:end-1,:)\x(2:end,:);
    evsA = eig(rhomat);
    rho(b)=sum(abs(evsA-1))/2;
end

[xs,ys]=kdfft1(rho,'knorm',1024,0.001);
figure,hold on
plot(xg*T,yg);
plot(xs*T,ys,'black')

prho = prctile(rho*T,95)
pdi = prctile(dis*T,95)
plot([1,1]*prho,[0,max(ys)])
plot([1,1]*pdi,[0,max(ys)])

sum(dis*T>prho)/M
sum(dis*T>pdi)/M

sum(rho*T>pdi)/B
sum(rho*T>prho)/B