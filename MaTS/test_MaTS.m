% skript for testing the alternating optimization approach for the
% estimation of autoregressions for matrix valued time series. 

T = 5000;
M = 3;
N = 2;
L = 2;
p = 2;

% exogenous part 
Mz = 4;
Nz = 4;
Lz = 2;
pz = 3; 

Sigma = eye(M*N);

% check normalize generated system
% test normalization 
% define matrices 
Co = zeros(M,Mz,L,p);
Co=rand(M,Mz,L,p);

Do = zeros(N,Nz,L,p);
Do=rand(N,Nz,L,p);

% calculate corresponding Kronecker product in vectorization 
vec_poly_before =  vectorized_syst(Co,Do);

% normalize 
[C,D] = norm_MaTS_syst(Co,Do);
vec_poly_after =  vectorized_syst(C,D);

% compare output: should be numerically zero. 
norm(vec_poly_before- vec_poly_after,"fro")


% specify system 

% estimate R systems to check, if it works 
R=100;

f = waitbar(0,'Please wait...');
est_err = zeros(R,4); 

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

    % exogenous part 
    if (Mz>0)
        C = randn(M,Mz,Lz,pz);
        D = randn(N,Nz,Lz,pz);
        [C,D] = norm_MaTS_syst(C,D);
    else
        C = [];
        D= [];
    end

    % generate data
    [Y,Z] = simu_MaTS(A,B,C,D,T,Sigma);

    % detect instabilities left
    while (norm(Y,"fro")>100000000)||(sum(sum(squeeze(sum(isnan(Y)))>0)))
        B = B/2;
        [A,B] = norm_MaTS_syst(A,B);
        [Y,Z] = simu_MaTS(A,B,C,D,T,Sigma);
    end

    vec_poly =  vectorized_syst(A,B);
    vec_poly_z =  vectorized_syst(C,D);
    % estimate 
    [Aest,Best,Cest,Dest]= est_MaTS(Y,L,p,Z,Lz,pz);

    % calculate vectorized system 
    vec_poly_est =  vectorized_syst(Aest,Best);
    vec_poly_est_z =  vectorized_syst(Cest,Dest);
    % compare output 
    est_err(r,1) = norm(A-Aest,"fro");
    est_err(r,2) = norm(B-Best,"fro");
    est_err(r,3) = norm(vec_poly-vec_poly_est,"fro");
    est_err(r,4) = norm(vec_poly_z-vec_poly_est_z,"fro");

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

    if (est_err(r,3)+est_err(r,4)>1)
        break
    end
end
close(f)

% see, if errors are small
hist(est_err(:,3:4))

