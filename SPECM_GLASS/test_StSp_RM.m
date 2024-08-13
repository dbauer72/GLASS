% test estimation using state space for each region. 

% specs 
N = 4; % number regions 
M = 3; % number of vars 
T = 2000; % number of time point 
ni=3;
n_full = N*ni; % number of states total (4 per region)
s_full= N*M; 

% define W matrices: 
% dominant variable average over all variables
W_w = randn(0,s_full);

% star variables: average over all vars in all other regions, specific to
% vars. 
Wst = kron(ones(1,N),eye(M))/3;

for n=1:N
    Wstar{n} = Wst;
    Wstar{n}(:,(n-1)*M+[1:M]) = 0; 
end

% draw system without spillovers 
A = zeros(n_full,n_full);
C = zeros(s_full,n_full);
K = C'; 
D = C; 
D_full = D; 
for n=1:N
    Ai= eye(ni)*0.5;

    while max(abs(eig(Ai)))>0.95
        Ai = randn(ni,ni);
    end
    Ci = eye(ni); %randn(M,ni); 
    Ki = eye(ni); %randn(ni,M);

    barAi= Ai - Ki*Ci;
    while max(abs(eig(barAi)))>0.95
        Ki = Ki/2;
        barAi= Ai - Ki*Ci;
    end

    ind = [(n-1)*ni+1:n*ni];
    indm = [(n-1)*M+1:n*M];
    A(ind,ind) = Ai;
    C(indm,ind) = Ci;
    K(ind,indm)=Ki; 
    D(indm,indm) = randn(ni,ni); 
end

Kostar = eye(n_full,s_full);
Costar = zeros(s_full,n_full);

for n=1:N
    Costar((n-1)*ni+[1:ni],(n-1)*ni+[1:ni]) = randn(ni,ni);
end

W = zeros(s_full,s_full);
for n=1:N
    W((n-1)*ni+[1:ni],:) = Wstar{n};
end

mv = max(abs(eig(A-Kostar*(Costar-W*C))));
while (mv>1)
    Kostar = Kostar/2;
    mv = max(abs(eig(A-Kostar*(Costar-W*C))));
end

mv = max(abs(eig(A-K*C-Kostar*Costar)))
while (mv>1)
    Kostar = Kostar/2;
    mv = max(abs(eig(A-K*C-Kostar*Costar)));
end

D_full = eye(s_full)-D*W;


% fill in system for dgp. 
A = A-Kostar*(Costar - W*C); 
K = K+Kostar*W*C; 

th = ss2ech_n(A,K,C); 

% generate observations 
u = randn(s_full,T);
x_full = zeros(n_full,1);
y_full = zeros(s_full,T);
for t=1:T
    y_full(:,t) = C*x_full + u(:,t);
    x_full = A*x_full + K*u(:,t);
end

y_full = D_full*y_full;

% cycle over regions to estimate systems
%wt = W_w*y_full;

% normalize to avoid issues with fit. 
%vv = wt*wt'/T;
%W_w = W_w/sqrt(vv);
%wt = wt/sqrt(vv);

indices = zeros(N,3);
ni = 6; 

for n=1:N
    yi = y_full([(n-1)*M+1:n*M],:);
    yistar = Wstar{n}*y_full;

    yif = [yi;yistar];
    indices(n,:) = [M,0,0];
    [resulti,thci,Aci,Kci,Cci,Omegaci,thii,thii2,llei] = SPECM_RM(yif',ni,indices(n,:),15,0);
    ths(n) = thci;
    results(n)=resulti; 
end

[th_stacked,ve_full,x_full] = stack_regional_models(ths,indices,Wstar,W_w,y_full,0);

