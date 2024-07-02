function [Aest,Best,Piest,Cest,Dest]= est_MaTS_OLS(Y,L,p,Z,Lz,pz)
% est_MaTS_OLS estimates the array pair (A,B) in a multi-component matrix valued AR(p)
% system. Estimation is performed via OLS in the vectorized system, then approximating the vectorized matrices. 
%
% SYNTAX:  [Aest,Best,Cest,Dest]= est_MaTS_OLS(Y,L,p,Z,Lz,pz,tol,maxit,printlevel)
%
% INPUT:   Y           ... T x M x N array of observations.
%          L           ... integer; lag length.
%          p           ... integer; number of terms per lag.
%          Z           ... T x Mz x Nz array of observations of exogenous vars.
%          Lz           ... integer; lag length for exo part.
%          pz           ... integer; number of terms per lag for exo part.
%
% OUTPUT:  Aest ... four dimensional array M x N x L x p. 
%          Best ... four dimensional array M x N x L x p.
%          Cest ... four dimensional array M x Mz x Lz x pz.
%          Dest ... four dimensional array Nz x N x Lz x pz.
% 
% AUTHOR: dbauer, 28.6.2024. 


if nargin<6
    pz = 0;
end
if nargin<5
    Lz = 0;
end
if nargin<4
    Z = [];    
end

Leff = max(L,Lz); 

dims = size(Y);
T = dims(1)-Leff; % effective sample size
M = dims(2);
N = dims(3); 
 
if isempty(Z)
    exo = 0;
    Lz = 0;
    pz =0;
    Mz = 0;
    Nz = 0;
else
    exo = 1;
    dimsz = size(Z);
    Mz = dimsz(2);
    if length(dimsz)>2
        Nz = dimsz(3);
    else
        Nz = 1;
    end
end

% vectorizations: predefine for numerical speed. 
% Warning: These matrices might get really large!
vY = zeros(T,M*N);
vYL = zeros(T,M*N*L);
if exo
    vZLz= zeros(T,Mz*Nz*(Lz+1));
else
    vZLz = zeros(T,0);
end

% fill in data in vectorized matrices 
for t=1:T
    vY(t,:)= reshape(Y(t+Leff,:,:),1,M*N);
    for j=1:L
        vYL(t,(j-1)*M*N+[1:M*N])= reshape(Y(t+Leff-j,:,:),1,M*N);
    end
    if exo
        vZLz(t,1:Mz*Nz) = reshape(Z(t+Leff,:,:),1,Mz*Nz);
        for jz=1:Lz
             vZLz(t,j*Mz*Nz+(1:Mz*Nz)) = reshape(Z(t+Leff-jz,:,:),1,Mz*Nz);
        end
    end
end

% calculate OLS 
OLS = ([vYL,vZLz]\vY)';

Aest = zeros(M,M,L,p);
Best = zeros(N,N,L,p);
Cest = zeros(M,Mz,Lz+1,pz);
Dest = zeros(N,Nz,Lz+1,pz);
Piest = zeros(M*N,Mz*Nz,Lz+1); % vectorized estimate
% obtain output 
for j=1:L
    vAj = OLS(:,(j-1)*N*M+(1:(N*M)));
    [Aj,Bj] = approx_terms(vAj,M,M,N,N,p);
    Aest(:,:,j,:) = Aj;
    Best(:,:,j,:)= Bj; 
end
if exo
    for jz=0:Lz   
        vAj = OLS(:,L*N*M+(jz)*Nz*Mz+(1:(Nz*Mz)));
        Piest(:,:,jz+1)= vAj;
        [Aj,Bj] = approx_terms(vAj,M,Mz,N,Nz,pz);
        Cest(:,:,jz+1,:) = Aj;
        Dest(:,:,jz+1,:)= Bj; 
    end
end


