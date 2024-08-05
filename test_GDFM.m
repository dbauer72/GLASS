% script to test functions 
% main integers: 
T = 1000;
N = 20; 
q = 3;
r = 6;
n=5;

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
Lambda = par2ortho_plt(randn(N*r),N,r)*sqrt(N);
th = theta_urs();
th.C = rand(r,n);
th.A = diag(0.8*(rand(n,1)-.5));
th.B = ones(n,q)+ randn(n,q)*0.1; 
th.D = rand(r,q); 
th.Omega = eye(q);

th_chi= th;
th_chi.C = Lambda*th.C;
th_chi.D = Lambda*th.D; 

% simulate the process
[y,chi,u]= simu_GDFM(T,ths,th,Lambda);

Gam = chi*chi'/T;
svd(Gam)
%%%%%%%%%%%%%%%%%%%%%%%%
% static PCA over time. 
%%%%%%%%%%%%%%%%%%%%%%%%
% 
GammaT = y*y'/T; 
[U,S] = svd(GammaT);
Ce_stat = U(:,1:(r+q)); 

plot(log(diag(S)),'x');

%%%%%%%%%%%%%%%%%%
% dynamic PCA
%%%%%%%%%%%%%%%%%%

% estimate spectrum pointwise. 
F =100;
freq = 0:(1/F):0.5;
BT = F; 

Sig_est = spec_mat(y,freq,BT);

% extract square root of dominant q cols.
LTs = rr_square_root(Sig_est,q);


% condense estimate into C matrix.
Ce = condense_square_root(LTs); 

% same for chi 
Sig_est_chi = spec_mat(chi,freq,BT);
LTs_chi = rr_square_root(Sig_est_chi,q);
Ce_chi = condense_square_root(LTs_chi); 

r=5;
k=4;
the2 = tf_est_spec(LTs_chi,freq,q,r,k);
the2 = stsp2inpnor(the2);

Ctr = [th_chi.C,th_chi.D];

LTs_F = zeros(size(Ctr,2),size(LTs_chi,2),BT/2+1);
for j=1:size(LTs_chi,3)
    LTs_F(:,:,j) = Ctr\squeeze(LTs_chi(:,:,j));
end
[thed,GammaTF] = tf_est_spec(LTs_F,freq,q,r,k);
thed.C = Ctr * thed.C;
thed.D = Ctr * thed.D; 
thed = stsp2inpnor(thed);

F = Ctr\chi;
Sig_est_F = spec_mat(F,freq,BT);
LTs_F2 = rr_square_root(Sig_est_F,q);
[thed2,GammaTF2] = tf_est_spec(LTs_F2,freq,q,r,k);
thed2.C = Ctr * thed2.C;
thed2.D = Ctr * thed2.D; 
thed2 = stsp2inpnor(thed2);





%%% try direct estimation with estimated covariance sequence
kmax = 20;

GammaT_est = zeros(size(F,1),size(F,1),kmax);
GammaT_est(:,:,1) = F*F'/T;

for jk=1:(kmax-1)
    GammaT_est(:,:,jk+1) = F(:,(jk+1):T)*F(:,1:(T-jk))'/T;
end;

[th_est_F] = CCA_sing_cov(GammaT_est,r,q,k,k,1);
th_est_F.C = Ctr * th_est_F.C;
th_est_F.D = Ctr * th_est_F.D; 
th_est_F = stsp2inpnor(th_est_F);

[the,Ae,Ke,CCe,De,Omegae] = CCA_sing(F',r,1,k,k,1);
the.C = Ctr * the.C;
the.D = Ctr * the.D; 
the = stsp2inpnor(the);

% direct estimation using CVA
[theX] = CCA_sing(y',r,1,4,4,1);
theX = stsp2inpnor(theX);

%%% estimation using F via spectrum
Fe = Ctr\y;
Sig_est_F = spec_mat(Fe,freq,BT);
LTs_F = rr_square_root(Sig_est_F,q,1);

r=5;
k=4;
theF = tf_est_spec(LTs_F,freq,q,r,k);
theF.C = Ctr * theF.C;
theF.D = Ctr * theF.D; 
theF = stsp2inpnor(theF);

%%%%%%%%%%%%%%%%%%
%% evaluate fit 
%%%%%%%%%%%%%%%%%%

% true TF
th_chi = stsp2inpnor(th_chi); 

% estimates: 
% the2     spectral for chi
% theF     spectral for F
% thed     reduced spectral for chi
% th_est_F cov based for F
% the      CVA for F
% theX     CVA for chi


% Omegas 
[th_chi.Omega,the2.Omega,theF.Omega,thed.Omega,th_est_F.Omega,the.Omega,theX.Omega]

% relative 
[th_chi.Omega,the2.Omega,theF.Omega,thed.Omega,th_est_F.Omega,the.Omega,theX.Omega]/th_chi.Omega

Ball = [th_chi.B,the2.B,theF.B,thed.B,th_est_F.B,the.B,theX.B]

Dall = [th_chi.D,the2.D,theF.D,thed.D,th_est_F.D,the.D,theX.D];
Dall(1:5,:)

Aall = [th_chi.A,the2.A,theF.A,thed.A,th_est_F.A,the.A,theX.A];

eAall = zeros(r,7); 
for j=1:7
    Ac = Aall(:,(j-1)*r+[1:r]);
    ev = eig(Ac);
    eAall(:,j)=ev(:);
end

eAall

Call = [th_chi.C,the2.C,theF.C,thed.C,th_est_F.C,the.C,theX.C];

Call(1:5,:)

%% impulse responses.

J=20; 
IM_diff = zeros(J+1,7);
for j=0:J
    IMall = [th_chi.C*th_chi.A^j*th_chi.B,the2.C*the2.A^j*the2.B,theF.C*theF.A^j*theF.B,thed.C*thed.A^j*thed.B,th_est_F.C*th_est_F.A^j*th_est_F.B,the.C*the.A^j*the.B,theX.C*theX.A^j*theX.B];
    IM_diff(j+1,:) = sum(abs(IMall - IMall(:,1)*[0,ones(1,6)]));
end
IM_diff


plot(IMall')

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate factors 
%%%%%%%%%%%%%%%%%%%%%%%%%%
Fest = Ce'*y; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% state space model for the factors 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[the,Ae,Ke,CCe,De,Omegae] = CCA_sing(Fest',5,1,5,5,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% combine to get model for common part chi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
the_chi = th_chi; 
the_chi.C = Ce*the.C; 
the_chi.A = the.A;
the_chi.B = the.K(:,1:q);
the_chi.D = Ce*De;
the_chi.Omega = De'*the.Omega*De;

% normalise result. 
the_chi2 = normalise_singular_stsp(the_chi); 

th_chi2 = normalise_singular_stsp(th_chi); 

% compare impulse response sequence. 
I = 10;
im_tr = zeros(n,q*I);
im_est = zeros(n,q*I);

im_tr(:,1:q)= th_chi2.D; 
im_est(:,1:q)= the_chi2.D; 

Cthtr = th_chi2.C; 
Cthe = the_chi2.C; 

Athtr = th_chi2.A; 
Athe = the_chi2.A;

Bthtr = th_chi2.B;
Bthe = the_chi2.B; 


for ij = 1:I
    im_tr(:,ij*q+[1:q]) = Cthtr * Bthtr;
    Bthtr = Athtr * Bthtr; 
    im_est(:,ij*q+[1:q]) = Cthe * Bthe;
    Bthe = Athe * Bthe; 
end

figure
mesh(im_tr-im_est)
