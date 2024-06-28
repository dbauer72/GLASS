% skript for testing the qMLE optimization approach for the
% estimation of autoregressions for matrix valued time series. 

T = 1000;
% M variables in N countries
M = 3;
N = 5;
% L lags in the AR  equation with p terms each 
L = 2;
p = 2;
% rr cointegrating relations, with Jr terms 
rr = 2;
Jr = 1; 
r = [rr,Jr]; 

% variance matrix with correlation in the cross section. 
Sigma = randn(M*N,M*N) + 5* eye(M*N); 

% estimate R systems to check, if it works 
R=10;

f = waitbar(0,'Please wait...');
est_err = zeros(R,4); 

% generate system 
A = zeros(M,M,L,p); % dimensions of A are M x N x max lag x number of components. 
A(:,:,1:L,1:p) = randn(M,M,L,p)*0.5;

B = zeros(N,N,L,p);
B(:,:,1:L,1:p) = randn(N,N,L,p)*0.5;

% stabilize to avoid explosive systems 
KroBA = zeros(M*N,M*N);
for j=1:L
    for k=1:p
       KroBA = KroBA + kron(squeeze(B(:,:,j,k))',squeeze(A(:,:,j,k)));
    end
end
lambda = max(abs(eig(KroBA)));
if (lambda>0.99)
    A = A/lambda * rand(1)*0.9;
end
[A,B] = norm_MaTS_syst(A,B);

% Pi matrix = alpha * beta'
[~,~,alpha] = param_term_rect(randn(100),M,1,N,rr,Jr,true);
[~,~,beta] = param_term_rect(randn(100),M,1,N,rr,Jr,false);
[alpha,beta] = norm_MaTS_VECM(alpha,beta,M,N,r)

% generate data
[Y] = simu_MaTS_VECM(A,B,alpha,beta,T,Sigma);

% detect instabilities left
while (norm(Y,"fro")>100000000)||(sum(sum(squeeze(sum(isnan(Y)))>0)))
    B = B/2;
    [A,B] = norm_MaTS_syst(A,B);
    beta = beta/2;
    [alpha,beta] = norm_MaTS_VECM(alpha,beta,M,N,r)
    [Y] = simu_MaTS_VECM(A,B,alpha,beta,T,Sigma);
end

% calculate true param for comparison
theta0 = MaTS2param_VECM(A,B,alpha,beta,r);
[A,B,alpha,beta] = param2MaTS_VECM(theta0,M,N,r,L,p)
theta0 = MaTS2param_VECM(A,B,alpha,beta,r);


vec_all = zeros(R,(N*M)*L);
th_all = zeros(R,length(theta0));

for jr=1:R 
    waitbar(jr/R,f,sprintf('Processing %d/%d',jr,R));
    [Y] = simu_MaTS_VECM(A,B,alpha,beta,T,Sigma);

    % estimate using qMLE
    [thetaest,Aest,Best,alphaest,betaest,ve]= qMLE_MaTS_VECM(Y,r,L,p);

    % compare output 
    [vec_poly] = vectorized_syst(A,B);    
    [vec_poly_est] = vectorized_syst(Aest,Best);
    for l=1:L
        vec_all(jr,(N*M)^2*(l-1)+[1:(N*M)^2])= reshape(squeeze(vec_poly_est(:,:,l)),1,(N*M)^2);
    end

    th_all(jr,:)= thetaest(:)';

    est_err(jr,1) = norm(A-Aest,"fro");
    est_err(jr,2) = norm(B-Best,"fro");
    est_err(jr,3) = norm(theta0-thetaest,"fro");
    est_err(jr,4) = norm(vec_poly-vec_poly_est,'fro');
    est_err(jr,5) = norm(alphaest*betaest'-alpha*beta','fro');

    for j=1:L
        Bmat = B(:,:,j,:);
        Bmatest = Best(:,:,j,:);
        if norm(Bmat-Bmatest,"fro")>norm(Bmat+Bmatest,"fro")
            Best(:,:,j,:) = -Best(:,:,j,:);
            Aest(:,:,j,:) = -Aest(:,:,j,:); 
        end
    end
    est_err(jr,1) = norm(A-Aest,"fro");
    est_err(jr,2) = norm(B-Best,"fro");

    if (est_err(jr,3)+est_err(jr,4)>100)
    %    break
    end
end
close(f)

% see, if errors are small
hist(est_err(:,4))

% compare variances to empirical variances
