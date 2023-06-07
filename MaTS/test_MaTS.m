% skript for testing the alternating optimization approach for the
% estimation of autoregressions for matrix valued time series. 

T = 10000;
M = 2;
N = 3;
L = 2;
p = 2;

Sigma = eye(M*N);

% check normalize generated system
% test normalization 
% define matrices 
Ao = zeros(M,M,L,p);
Ao=rand(M,M,L,p);

Bo = zeros(N,N,L,p);
Bo=rand(N,N,L,p);

% calculate corresponding Kronecker product in vectorization 
vec_poly_before =  vectorized_syst(Ao,Bo);

% normalize 
[A,B] = norm_MaTS_syst(Ao,Bo);
vec_poly_after =  vectorized_syst(A,B);

% compare output: should be numerically zero. 
norm(vec_poly_before- vec_poly_after,"fro")


% specify system 

% estimate R systems to check, if it works 
R=500;

f = waitbar(0,'Please wait...');
est_err = zeros(R,3); 

for r=1:R 
    waitbar(r/R,f,sprintf('Processing %d/%d',r,R));
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
    Y = simu_MaTS(A,B,T,Sigma);

    % detect instabilities left
    while (norm(Y,"fro")>100000000)||(sum(sum(squeeze(sum(isnan(Y)))>0)))
        B = B/2;
        [A,B] = norm_MaTS_syst(A,B);
        Y = simu_MaTS(A,B,T,Sigma);
    end

    vec_poly =  vectorized_syst(A,B);

    % estimate 
    [Aest,Best]= est_MaTS(Y,L,p);

    % calculate vectorized system 
    vec_poly_est =  vectorized_syst(Aest,Best);

    % compare output 
    est_err(r,1) = norm(A-Aest,"fro");
    est_err(r,2) = norm(B-Best,"fro");
    est_err(r,3) = norm(vec_poly-vec_poly_est,"fro");

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

    if (est_err(r,3)>1)
        break
    end
end
close(f)

% see, if errors are small
hist(est_err(:,3))
