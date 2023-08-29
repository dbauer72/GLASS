function [crit,grad,Hess,JV,ve,psi] = cal_crit_MaTS(Y,theta,L,p)
% cal_crit calculates the concentrated (w.r.t. the parameters for the variance) 
% Gaussian likelihood criterion function.
%
% SYNTAX: [crit,ve,psi] = cal_crit_MaTS(Y,theta,M,N,L,p)
%
% INPUTS: Y       ... T x M x N real matrix of observations.
%         theta   ... d vector of parameters
%         M       ... integer; 
%         N       ... integer
%         L       ... integer; number of lags
%         p       ... integer; number of terms per lag.
%
% OUTPUTS: crit   ... real, value of criterion function.
%          ve     ... M*N x T real matrix of vectorized residuals.
%          psi    ... M*N x T x d real matrix of derivative of residuals.
%
% AUTHOR: dbauer, 21.8.2023.

% get dimensions.
[T,M,N] = size(Y); 
nth = length(theta);
grad= zeros(nth,1);
Hess = zeros(nth,nth); 

% calculate AR matrices from parameters 
if nargout>1
    [A,B,dA,dB] = param2MaTS(theta,M,N,L,p);
else
    [A,B] = param2MaTS(theta,M,N,L,p);
end

% set up vectorized system 
% for each lag a M*N x M*N matrix. 
vPsi = zeros(M*N,M*N,L); 
for l=1:L
    for jr=1:p
        vPsi(:,:,l) = vPsi(:,:,l) + kron(squeeze(B(:,:,l,jr)),squeeze(A(:,:,l,jr)));
    end
end

if nargout >1
    dvPsi = zeros(M*N,M*N,L,length(theta)); 
    for d=1:length(theta)
        for l=1:L
            for jr=1:p
                dvPsi(:,:,l,d) = dvPsi(:,:,l,d) + kron(squeeze(B(:,:,l,jr)),squeeze(dA(:,:,l,jr,d))) + kron(squeeze(dB(:,:,l,jr,d)),squeeze(A(:,:,l,jr)));
            end
        end
    end
end


% calculate residuals 
ve = zeros(T,M*N); 
% initial vallues for residuals set to NaN
ve(1:L,:) = NaN;  
%Omega = zeros(M*N,M*N); 
%for t=(L+1):T % cycle over time 
%    vYt = reshape(squeeze(Y(t,:,:)),M*N,1);
%    vet = vYt; 
%    for l=1:L % cycle over 
%        Ytml = reshape(squeeze(Y(t-l,:,:)),M*N,1);
%        vet= vet - squeeze(vPsi(:,:,l)) * Ytml;
%    end
%    Omega = Omega + vet*vet'; 
%    ve(t,:)= vet; %reshape(vet,M,N);
%end

%Omegahat = Omega/(T-L); 

% alternative: vectorize
vY = reshape(Y,T,M*N);

ve2 = vY((L+1):T,:);

for l=1:L
    ve2 = ve2 - vY((L+1-l):(T-l),:)*squeeze(vPsi(:,:,l))';
end

Omegahat = ve2'*ve2/(T-L);

% calculate criterion value 
crit = log(det(Omegahat));

% calculate derivative of residuals if needed
if nargout>1
    psi = zeros(T,N*M,nth);
    iOm = inv(Omegahat);
    for d=1:length(theta)
        %for t=(L+1):T % cycle over time
        %    vYt = reshape(squeeze(Y(t,:,:)),M*N,1);
        %    vet = vYt;
        %    psit = vet*0;
        %    for l=1:L % cycle over
        %        Ytml = reshape(squeeze(Y(t-l,:,:)),M*N,1);
        %        psit= psit - squeeze(dvPsi(:,:,l,d)) * Ytml;
        %    end
        %    psi(t,:,d)=psit;
        %end

        % alternative: vectorized version
        psi2 = ve2*0;
        for l=1:L
             psi2= psi2 - vY((L+1-l):(T-l),:)*squeeze(dvPsi(:,:,l,d))';            
        end
        % gradient calc: grad_i = tr[ Omega^{-1} dOmega]
        %grad(d) = 2*trace(iOm * ve((L+1):end,:)' * squeeze(psi((L+1):end,:,d)))/(T-L);
        grad(d) = 2*trace(iOm * ve2' * psi2)/(T-L); 
        psi((L+1):end,:,d)=psi2; 
    end

    if nargout>2
        % Hess cycle over twice times d. 
        % same for JV.
        for d=1:length(theta)
            psid = squeeze(psi((L+1):end,:,d));
            dOmegad = ve2'*psid /(T-L);
            for e=1:d
                % second derivative of ve 
                d2psi = ve2*0;
                for l=1:L
                    d2vPsi = zeros(M*N,M*N); 
                    for jr=1:p
                        d2vPsi = d2vPsi + kron(squeeze(dB(:,:,l,jr,e)),squeeze(dA(:,:,l,jr,d))) + kron(squeeze(dB(:,:,l,jr,d)),squeeze(dA(:,:,l,jr,e)));
                    end
                    d2psi= d2psi - vY((L+1-l):(T-l),:)*d2vPsi';            
                end
    
                psie = squeeze(psi((L+1):end,:,e));
                dOmegae = ve2'*psie /(T-L);
                Hess(d,e) = 2* trace( iOm * (psie' * psid + d2psi'*ve2) )/(T-L)-2*trace( iOm * dOmegad * iOm * dOmegae);
                JV(d,e) = 4* trace( iOm * (psie' * psid) )/(T-L)^2;
            end
        end   

        % take symmetry into account 
        dH = diag(Hess);
        Hess = Hess + Hess' - diag(dH);

        dJ = diag(JV);
        JV = JV + JV' - diag(dJ);
end
    
end

