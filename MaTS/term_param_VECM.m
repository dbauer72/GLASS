function [theta] = term_param_VECM(alpha,M,N,p,which);
% term_param calculates the parameters corresponding to the p terms in (A,B).
%
% SYNTAX: [theta] = term_param_VECM(A,B,r,which);
%
% INPUTS: alpha ... M*N x r real array of p terms.
%         M,N ... integers; dimensions of observations. 
%         Jr  ... integers; number of terms.
%         which ... Boolean; true: parameters for alpha. else: params for
%                   beta. 
%        
% OUTPUTS: theta ... d vector of parameters.
%
% REMARKS: used in VECM version of MaTS Y_t = \sum_{j=1}^p A_j Y_{t-1} B_j' + E_t.
% provides parameters (using the QR decomposition) 
% for r(2) terms matrices.
% 
% AUTHOR: dbauer, 21.8.2023.

% dimensions 
[MN,rr] = size(alpha); 
if (MN ~= M*N)
    error('term_param_VECM: dimensions do not match!');
end

% transform using Phi.
Phi = reorder_KronProd(alpha,N,M,rr,1);

% calculate QR decomposition 
[Q,R]= qr(Phi');

% take out the relevant matrices
Q = Q(:,1:p);
R = R(1:p,:);


% adjust sign of diagonal entries
signs = sign(diag(R(1:p,1:p)));
Q(:,1:p)=Q(:,1:p)*diag(signs);
R(1:p,:)=diag(signs)*R(1:p,:);

% obtain parameters
theta = ortho2par_LR(Q);

% parameters for the R part. 
R_A = NaN(N*p,rr);
for j=1:p
   R_A((j-1)*N+[1:N],:)=reshape(R(j,:),rr,N)';
end

% reorder to impose normalisation assumptions  
ind = 1:(N*p); 
N_ind = (reshape(ind,N,p))';
R_A = R_A(N_ind(:),:);

% impose assumptions by introducing NaNs
for j=2:p
    R_A(j,1:(j-1))=NaN;
end

if which
    R_A(end-rr+[1:rr],:)=NaN;
end

% vectorize 
vR_A = R_A(:);
% find non-NaNs 
ind = find(isnan(vR_A)==0);
theta_A = vR_A(ind);

% combine 
theta = [theta(:);theta_A(:)];

% if which % parameters for alpha: first row starts with vectorization of identity matrix. 
%     lIr = rr^2; 
%     if lIr<p
%         theta(end+[1:(p-lIr)])=R(1,(lIr+1):p);
%     else
%         lIr = lIr - p;
%     end
% 
%     for jr=2:p
%         theta(end+[1:(p-jr+1)])=R(jr,jr:p);
%     end
% 
%     th = reshape(R(:,(p+1):end)',(M*rr-p)*p,1);
%     theta = [theta(:);th(lIr+1:end)];
% 
% else 
%     for jr=1:p
%         theta(end+[1:(p-jr+1)])=R(jr,jr:p);
%     end
%     theta = [theta(:);reshape(R(:,(p+1):end)',(M*rr-p)*p,1)];
% end
