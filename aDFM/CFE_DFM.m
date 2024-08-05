% sample size 
T = 200;
N = 200;
% number of repetitions 
M = 1000; 

% define system 
% idiosyncratic component: AR(1) with random rho. 
ths(1) = theta_urs();
ths(1).A = 1.4*(rand(1)-.5);
ths(1).B = 1;
ths(1).C = ths(1).A;
ths(1).D = 1;
ths(1).Omega = .5; 
for j=2:N
    ths(j) = theta_urs();
    ths(j).A = 1.4*(rand(1)-.5);
    ths(j).B = 1;
    ths(j).C = ths(j).A;
    ths(j).D=1;
    ths(j).Omega = .5; 
end

% common factor part. 
n = 3; 
r = 5; 
q = 2;
p=2;
th_chi = theta_urs();
th_chi.C = randn(r,n);
th_chi.A = diag(0.8*(rand(n,1)-.5));
th_chi.A = diag([-.8,.4,.8]);
%dA = diag(th_chi.A);
%dA(abs(dA)<0.5)=0.5;
%th_chi.A = diag(dA); 
th_chi.B = randn(n,q); 
th_chi.D = randn(r,q); 
th_chi.Omega = eye(q);

% calculate the transfer function for differenced series
dth_chi = differentiate_StSp_sing(th_chi);
ith_chi = integrate_StSp_sing(th_chi);

% loadings 
Lambda = randn(N,r); 


% start the engine. 
systs = zeros(M,n^2+n*r+n*q+r*q,3);
integers = zeros(M,5,3);
sing_vals = zeros(M,10,3);
norms = zeros(M,11,3);

% demonstrate eigenvalue property 
[y,chi,u]= simu_GDFM(4000,ths,th_chi,Lambda);

Gammaf = y*y'/T; 
evs = zeros(N,N);
for j=1:N
    evj = eig(Gammaf(1:j,1:j));
    evs(j,1:j)=sort(evj,'descend');
end

figure;
hold on;
for j=1:4
    plot(1:N,(evs(:,j)));
end;

set(gca,'fontsize',16)
xlabel('N')
ylabel('Eigenvalues')

print -dpng large_evs.png

figure;
hold on;
for j=5:200
    plot(j:N,(evs(j:N,j)),'.');
end;

set(gca,'fontsize',16)
xlabel('N')
ylabel('Eigenvalues')

print -dpng small_evs.png



me_no = zeros(4,4,11,3);
in_per = zeros(4,4,4,3);

me_phat = zeros(4,4,3); 

h = waitbar(0,"PLease wait ...");
Ts = [200,400,1600,800];
Ns = [50,100,150,200];
for nn=1:4
    N = Ns(nn);
    for tt= 1:4
        T = Ts(tt); 
        for m=1:M
            waitbar(m/M,h,sprintf("N: %d,T: %d, %d of %d",N,T,m,M));
            % generate data 
            % simulate the process
            [y,chi,u,x,e,F]= simu_GDFM(T,ths,th_chi,Lambda);
            
            % three different estimates: 
            % not overdifferenced 
            [th,rest,Lambdahat,pest,nest,nest2,s,qest,Omega] = cal_est_GDFM(y(1:N,:),r,n,q);
            if (length(s)>=2*n)
                sing_vals(m,1:2*n,1)=s(1:2*n); 
            else
                sing_vals(m,1:length(s),1)=s;
            end
            
            %systs(m,:)= [th.A(:)',th.B(:)',th.C(:)',th.D(:)'];
            HT = pinv(Lambda(1:N,:))*Lambdahat;
            norms(m,:,1)= cal_norm_dev(th,th_chi,HT,10); 
            integers(m,:,1)= [rest,pest,nest,nest2,qest];

            % overdifferenced 
            dy = y(1:N,2:end)-y(1:N,1:end-1); 
            [th,rest,Lambdahat,pest,nest,nest2,s,qest,Omega] = cal_est_GDFM(dy,r,n+q,q);
            if (length(s)>=2*n)
                sing_vals(m,1:2*n,2)=s(1:2*n); 
            else
                sing_vals(m,1:length(s),2)=s;
            end
            
            %systs(m,:)= [th.A(:)',th.B(:)',th.C(:)',th.D(:)'];
            HT = pinv(Lambda(1:N,:))*Lambdahat;
            norms(m,:,2)= cal_norm_dev(th,dth_chi,HT,10); 
            integers(m,:,2)= [rest,pest,nest,nest2,qest];

            % integrated 
            ichi = cumsum(chi');
            iy = ichi'+u; 
            
            [th,rest,Lambdahat,pest,nest,nest2,s,qest,Omega] = cal_est_GDFM(iy(1:N,:),r,n+q,q);
            if (length(s)>=2*n)
                sing_vals(m,1:2*n,3)=s(1:2*n); 
            else
                sing_vals(m,1:length(s),3)=s;
            end
                        
            HT = pinv(Lambda(1:N,:))*Lambdahat;
            dth = differentiate_StSp_sing(th);
            norms(m,:,3)= cal_norm_dev(dth,th_chi,HT,10); 
            integers(m,:,3)= [rest,pest,nest,nest2,qest];
        end
        for j=1:3
            me_no(nn,tt,:,j)= mean(squeeze(norms(:,:,j)));
            nj = 0; 
            switch j
                case 1
                    nj = n;
                case 2
                    nj = n+q;
                case 3
                    nj = n+q;
            end
            in_per(nn,tt,:,j)= [sum(integers(:,1,j)==r),sum(integers(:,2,j)==p),sum(integers(:,4,j)==nj),sum(integers(:,5,j)==q)]/M;
            me_phat(nn,tt,j) = mean(squeeze(integers(:,2,j)));
        end
    end
end

close(h); 
save CFE_simul_r5.mat

%%% plots for integers 
% rhat 
figure;
hold on;
hist(integers(:,1,:),[0:5])

set(gca,'fontsize',16)
xlabel('r')
ylabel('Frequency')
legend('stationary','overdifferenced','integrated')
title('Estimating r')

print -dpng rhat.png

% phat 
figure;
hold on;
hist(integers(:,2,:),[0:5])

set(gca,'fontsize',16)
xlabel('p')
ylabel('Frequency')
legend('stationary','overdifferenced','integrated')
title('Estimating p')

print -dpng phat.png

% mean p values 
me_phat(:,[1,2,4,3],:)


% nhat 
figure;
hold on;
hist(integers(:,4,:),[0:6])

set(gca,'fontsize',16)
xlabel('n')
ylabel('Frequency')
legend('stationary','overdifferenced','integrated')
title('Estimating n')

print -dpng nhat.png



% canonical corr. close to 1
figure;
hold on;
hist(integers(:,3,:),[0:4])

set(gca,'fontsize',16)
xlabel('n')
ylabel('Frequency')
legend('stationary','overdifferenced','integrated')
title('Can. Corr. Close to 1')

print -dpng ccahat.png

% percentages N=50, all four sample sizes, can. corr. close to 1. 
squeeze(in_per(1,[1,2,4,3],3,:))'


% qhat 
figure;
hold on;
hist(integers(:,5,:),[0:5])

set(gca,'fontsize',16)
xlabel('q')
ylabel('Frequency')
legend('stationary','overdifferenced','integrated')
title('Estimating q')

print -dpng qhat.png

% 
squeeze(in_per(:,[1,2,4,3],5,3))

% accuracy 
acc = squeeze(sum(me_no(1,[1,2,4,3],:,:),3));
figure;
hold on;
plot(Ts([1,2,4,3]),acc)

set(gca,'fontsize',16)
xlabel('T')
set(gca,'XTick',sort(Ts))
ylabel('Accuracy')
legend('stationary','overdifferenced','integrated')

print -dpng acc.png


% 
% % integers
% prct_correct = [sum(integers(:,1)==r),sum(integers(:,2)==p),sum(integers(:,3)==n),sum(integers(:,4)==q)]/M
% 
% % accuracy of system estimates 
% me = mean(systs);
% systs_dm = systs - ones(M,1)*me; 
% 
% mae = sum(abs(systs_dm))/M
% mse = sum(abs(systs_dm).^2)/M
% 
% % this is really bad. 
% k=3; 
% Hfp =  zeros(r*k,q*k);
% Hfpe = Hfp;
% 
% for a=1:k
%     for b=1:k
%         Hfp((a-1)*r+[1:r],(b-1)*q+[1:q]) = th_chi.C*th_chi.A^(a+b-1)*th_chi.B;
%         Hfpe((a-1)*r+[1:r],(b-1)*q+[1:q]) = th.C*th.A^(a+b-1)*th.B;
%     end
% end
% 
% % compare eigenvalues for A
% evals = zeros(M,5);
% for m=1:M
%     ev = sort(abs(eig(reshape(systs(m,1:25),5,5))),'ascend');
%     evals(m,:)= ev(:)';
% end
% 
% 


