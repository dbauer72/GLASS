% skript for testing the qMLE optimization approach for the
% estimation of autoregressions for matrix valued time series. 

T = 1000;
M = 3;
N = 2;
L = 2;
p = 2;

Sigma = randn(M*N,M*N) + 5* eye(M*N); 


% estimate R systems to check, if it works 
R=100;

f = waitbar(0,'Please wait...');
est_err = zeros(R,4); 

% only one system to compare variance estimation
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

    % generate data
    [Y] = simu_MaTS(A,B,[],[],T,Sigma);

    % detect instabilities left
    while (norm(Y,"fro")>100000000)||(sum(sum(squeeze(sum(isnan(Y)))>0)))
        B = B/2;
        [A,B] = norm_MaTS_syst(A,B);
        [Y] = simu_MaTS(A,B,[],[],T,Sigma);
    end

    % calculate true param for comparison
    theta0 = MaTS2param(A,B);

    vec_all = zeros(R,36*2);
    th_all = zeros(R,length(theta0));

for r=1:R 
    waitbar(r/R,f,sprintf('Processing %d/%d',r,R));
    [Y] = simu_MaTS(A,B,[],[],T,Sigma);

    % estimate using qMLE
    [thetaest,Aest,Best,Vest,sdA,sdB,sPsi,thetai,Ai,Bi]= qMLE_MaTS(Y,L,p);

    
    % compare output 
    [vec_poly] = vectorized_syst(A,B);    
    [vec_poly_est] = vectorized_syst(Aest,Best);
    for l=1:L
        vec_all(r,(N*M)^2*(l-1)+[1:(N*M)^2])= reshape(squeeze(vec_poly_est(:,:,l)),1,(N*M)^2);
    end

    th_all(r,:)= thetaest(:)';

    est_err(r,1) = norm(A-Aest,"fro");
    est_err(r,2) = norm(B-Best,"fro");
    est_err(r,3) = norm(theta0-thetaest,"fro");
    est_err(r,4) = norm(vec_poly-vec_poly_est,'fro');

    for j=1:L
        Bmat = B(:,:,j,:);
        Bmatest = Best(:,:,j,:);
        if norm(Bmat-Bmatest,"fro")>norm(Bmat+Bmatest,"fro")
            Best(:,:,j,:) = -Best(:,:,j,:);
            Aest(:,:,j,:) = -Aest(:,:,j,:); 
        end
    end
    est_err(r,1) = norm(A-Aest,"fro");
    est_err(r,2) = norm(B-Best,"fro");

    if (est_err(r,3)+est_err(r,3)>100)
        break
    end
end
close(f)

% see, if errors are small
hist(est_err(:,4))
