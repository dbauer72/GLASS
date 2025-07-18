% test BLT system estimation 

% integers:
index = [2,1];
sf = [3,2];
nf = [4,2];
n = nf(1); 

% system for star vars 
As = eye(2);
As(2,2)= 0.2*randn(1);
Cs = eye(2);

Ks = diag([0.5,0.4]);

% generate system
A = eye(nf(1));
A(3:4,3:4)=0.25*randn(2,2);

C = [eye(index(1)),randn(index(1),n-index(1));randn(sf(1)-index(1),n)];
K = randn(n,sf(1));

mA = abs(eig(A-K*C));
while (max(abs(mA))>0.99)
    K  = randn(n,sf(1));
    mA = abs(eig(A-K*C));
end

% add D and B 
D = randn(sf(1),sf(2));
B = randn(nf(1),sf(2));
B(:,1:index(2))=0; % integrated coordinate does not run through state.

LO= randn(sf(1),sf(1));
Omega = LO * LO';

% write into th structure.
th = theta_urs;
th.which = 'SS';
th.A = A;
th.K=K;
th.C=C;
th.D = D;
th.B = B;
th.Omega = Omega; 

% -- extract params 
[par,parI1] = th2param_BLT(th,index,1);

% --- convert param2mat
th2 = param2th_BLT(par,nf(1),sf,index);

par2 = th2param_BLT(th2,index,1);

th2 = param2th_BLT(par2,nf(1),sf,index);

% --- now generate signals ---
T = 1000;
% --- ystar: s=2,n=2,c=1. 
ystar = zeros(T,sf(2));
u = randn(T,sf(2));
ystar(:,1)=cumsum(u(:,1));
ystar(1,2)= u(1,2);
for t=2:T
    ystar(t,2)= ystar(t-1,2)*0.5+u(t,2);
end

plot(ystar)

% --- generate y ---
% th: filters y_t* and u_t.
y1 = idsim_x(T,th,ystar);

y = [ystar,y1];

% --- test likelihood calc ---
[qlike,tres] = cal_quasi_like_BLT(par,y,n,sf,index);
plot(tres)
Omh = tres' * tres/T;

% --- test RH approach ---
Abar = th.A-th.K*th.C

[LL,alphahat,betahat,Chat] = RH_specm_BLT(y,Abar,th.B-th.K*th.D,th.K,1,1,'n',5);

% --- test exogeneity testing 
[k,Omega,AICs,exo_test,tharx,thj,thst] = aicest_exo(y,3,20,0);
thi_AR = project_init_ARX(tharx,3,2,4,2,2,1);
parami = th2param_BLT(thi_AR,index,1);
options = optimoptions('fminunc','display','final');
options.MaxFunctionEvaluations = 10000;
%restrict.scale = ones(length(parami),1);
restrict.det_res = 0;

% estimate initial stable model
pare = est_cal_like_hess_BLT(y,4,sf,index,thi_AR,0,restrict);
resulte = compile_results_BLT(pare,4,sf,index,y,0,4);


% --- now test estimation ---
[result,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM_BLT(y,n,sf,index,5);

% --- include a model for the star vars
[resultst,thcst,Acst,Kcst,Ccst,Omegacst] = SPECM_I1(y(:,1:2),2,2,1,5,0);

%---- combine models ----
th_GL = combine_models_star(thc,thcst);

% --- compared to a full model for both 
[resultj,thj,Aj,Kj,Cj,Omegaj] = SPECM_I1(y,5,6,3,6,0);   
thj2 = project_init(thj,3,2,4);
restrict.det_res = 0;
parami = th2param_BLT(thj2,index,1);
options = optimoptions('fminunc','display','final');
options.MaxFunctionEvaluations = 10000;
%restrict.scale = ones(length(parami),1);

% estimate initial stable model
pare = est_cal_like_hess_BLT(y,4,sf,index,thj2,0,restrict);
resulte = compile_results_BLT(pare,4,sf,index,y,0,4);


%[parj2,Ak,Kk,Ck] = th2param(thj2,3,1);
%paro = extr_lowtri(thj2.Omega);
%restrict.det_res = 0;
%[Ae,Ke,Ce,De,Omegae] = param2syst([paro(:);parj2(:)],5,6,0,3,restrict);

%
%[qlike_GL,tres_GL] = cal_quasi_like([paro(:);par_GL(:)],y,5,0,6,3,0);


% ---- calculate parameters and evaluate likelihood
parj = th2param(thj,3,0);
paro = extr_lowtri(thj.Omega);

[qlikej,tresj] = cal_quasi_like([paro(:);parj(:)],y,5,0,6,3,0);

% --- same for combined system 
[par_GL,Ak,Kk,Ck] = th2param(th_GL,3,1);
paro_GL = extr_lowtri(th_GL.Omega);

[Ae,Ke,Ce,De,Omegae] = param2syst([paro(:);par_GL(:)],5,6,0,3,restrict);

restrict.det_res = 0;

[qlike_GL,tres_GL] = cal_quasi_like([paro(:);par_GL(:)],y,5,0,6,3,0);
% -- adjust Omega to enhance estimate
Omh = tres_GL'*tres_GL/T;
paro_GL = extr_lowtri(Omh);
[qlike_GL,tres_GL] = cal_quasi_like([paro_GL(:);par_GL(:)],y,5,0,6,3,0);

% --- now compare:
[qlikej,qlike_GL]

% --- number of parameters 
[length(parj),length(result.param)+length(resultst.param)]

% --- to true model 
th_comb = combine_models_star(th,thcst);

[par_comb] = th2param(th_comb,3,1);
paro_comb = extr_lowtri(th_comb.Omega);
restrict.det_res = 0;

[qlike_comb,tres_comb] = cal_quasi_like([paro_comb(:);par_comb(:)],y,5,0,6,3,0);

% --- compare AR representation for separate model to the one for the joint
% model:

% joint model 
Abar = th_GL.A- th_GL.K*th_GL.C; 
AR = zeros(3,5,20);
for j=1:20
    AR(:,:,j) = -[th_GL.C(3:5,:)-th.D*th_GL.C(1:2,:)]*Abar^(j-1)*th_GL.K;
end

% conditional model: 
Abarc = thc.A- thc.K*thc.C; 
ARc = zeros(3,5,20);
for j=1:20
    ARc(:,:,j) = -thc.C*Abarc^(j-1)*[thc.B-thc.K*thc.D,thc.K];
end

norm(AR-ARc,'fro')





