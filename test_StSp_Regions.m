% test estimation using state space for each region. 

% specs 
N = 5; % number regions 
M = 6; % number of vars 
T = 1000; % number of time point 
n = 6; % number of states total 


% draw system without star variables. 
A= randn(n,n);

while max(abs(eig(A)))>0.95
    A = randn(n,n);
end

C = randn(N*M,n); 
K = randn(n,N*M);

barA= A - K*C;
while max(abs(eig(barA)))>0.95
    K = K/2;
    barA= A - K*C;
end

th = ss2ech_n(A,K,C); 

% generate observations 
u = randn(N*M,T);
x = zeros(n,1);
y = zeros(N*M,T);
for t=1:T
    y(:,t) = C*x + u(:,t);
    x = A*x + K*u(:,t);
end

% cycle over regions to estimate systems
for j=1:N
    %ystar = Wi(j,:)*y;
    %ycom = Wom*y;

    %z = [y(j,:)',ystar'];
    z = y((j-1)*M+(1:M),:)';
    [k] = aicest(z,M,10);
    k = max(k,1);
    [thi(j)] = CCA(z,n,2*k,2*k,[]);
%    [result(j),thc,Ac,Kc,Cc,Omegac,thi,lle] = StSp_I0(z,M,n,n,1,thi(j));   
end


% combine states to large one
xh = zeros(T,n*N); 
res = u'*0;
for j=1:N
    % take from results
    %Ae = result(j).theta.A;
    %Ke = result(j).theta.K;
    %Ce = result(j).theta.C;
    % alternative: CCA estimates
    Ae = thi(j).A;
    Ke = thi(j).K;
    Ce = thi(j).C;
    
    Abar = Ae-Ke*Ce;
    xh(:,(j-1)*n+(1:n)) = ltitr(barA,[Ke],y((j-1)*M+(1:M),:)',zeros(n,1));
    res(:,(j-1)*M+(1:M)) = y((j-1)*M+(1:M),:)'-xh(:,(j-1)*n+(1:n))*Ce';
%    [res(:,(j-1)*M+(1:M)),~] = est_initial_val(res(:,(j-1)*M+(1:M)),barA,Ke,Ce);
end

% PCA on the state to eliminate not used components. 
V_x = xh'*xh/T;
[U,S,V] = svd(V_x);
plot(diag(S),'x')

% reestimate large system, based only on the state
x_red = xh*U(:,1:n);
Cest = (x_red\y')';
res_full = y' - x_red*Cest';
Aest = (x_red(1:end-1,:)\x_red(2:end,:))';
Kest = (res_full(1:end-1,:)\x_red(2:end,:))';

thest = ss2ech_n(Aest,Kest,Cest);

% compare results
ML = 10; 
IM = impulse(th,ML);
IMest = impulse(thest,ML);
for j=1:ML-1
    diff_norm_IM(j)=norm(IM(:,:,j+1)-IMest(:,:,j+1));
end

norm(IM(:,:,2:end),'fro')
diff_norm_IM

norm(res_full-u')/norm(u)

est_diff = res_full'*res_full/T - u*u'/T;
reshape(diag(est_diff),M,N)

figure; 
hold on; 
for j=1:N*M
    plot(res_full(:,j),u(j,:),'.');
end
