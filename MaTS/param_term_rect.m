function [A,B,alpha] = param_term_rect(theta,M1,M2,N1,N2,p,which)
% param_term converts the parameter vector into the sum of r terms in a
% matrix time series model
%
% SYNTAX:  [A,B,da,dB] = param_term(theta,M1,M2,N1,N2,p,which)
%
% INPUTS:  theta ... d dimensional real vector of parameters.
%          M1,M2     ... integers, dimensions for row matrices.
%          N1,N2     ... integers, column matrix dimensions.
%          p     ... integer, number of terms.
%          which ... Boolean. true indicates that parameters for alpha are
%                    supplied subject to identifying restrictions. 
%
% OUTPUTS: A  ... M1 x M2 x p array of real matrices. 
%          B  ... N1 x N2 x p array of real matrices.
%
% REMARKS: 
% for alpha uses a normalisation based on setting the tail matrix of stacked row matrices to identity.
% This works for (M-1)*Jr >= rr, in which case it is not inflicting with restrictions due to the QR factorization.  
% 
% AUTHOR: dbauer, 21.8.2023.

nth = length(theta);
Msq = M1*M2; 
ntriang = p*(p-1)/2;
npar = p*(Msq-p)+ntriang;
par = theta(1:npar); 
theta(1:npar)=[];
% calculate Q.
[Q] = par2ortho_LR(par,Msq,p);
for jr = 1:p
    A(:,:,jr) = reshape(Q(:,jr),M1,M2);
end

% calculate R
% if which is true, the R matrix starts with an identity. 
if which 
    Ir = eye(p);
    vIr = Ir(:);
end

% now fill in the R part.
R_A = NaN(N1*p,N2);
if which 
    R_A(end-N2+[1:N2],:)=eye(N2);
end
if p>1
    for j=2:N2
        R_A(j,1:(j-1))=0;
    end
end

% vectorize 
vR_A = R_A(:);
% find NaNs 
ind = find(isnan(vR_A));
vR_A(ind)= theta(1:length(ind));

% make matrix out of this 
R_A = reshape(vR_A,N1*p,N2);

% reorder to obtain alpha_{j,r} matrices. 
ind = 1:(N1*p); 
N_ind = (reshape(ind,p,N1))';
R_A = R_A(N_ind(:),:);

for jr=1:p
    B(:,:,jr) = R_A((jr-1)*N1+[1:N1],:)';
end


% R = zeros(p,N1*N2);
% for j=1:p
%     if (which && (~isempty(vIr)) && (j==1))
%         if length(vIr)<(p-j+1)
%             vIr((end+1):(p-j+1)) = theta(1: (p-j+1-length(vIr)));
%             theta(1: (p-j+1-length(vIr))) = [];
%         end
%         R(j,j:p) = vIr(1:(p-j+1));
%         vIr(1:(p-j+1)) = []; 
%     else
%         R(j,j:p) =theta(1:(p-j+1));
%         theta(1:(p-j+1)) = []; 
%     end
% end
% 
% if (which && (~isempty(vIr)))
%     theta = [vIr;theta(:)];
% end
% R(:,(p+1):end) = reshape(theta(1:((N1*N2-p)*p)),(N1*N2-p),p)';
% for jr = 1:p
%     B(:,:,jr) = reshape(R(jr,:)',N2,N1);
% end

alpha = zeros(M1*N1,M2*N2); 

for jr = 1:p
    alpha = alpha + kron(squeeze(B(:,:,jr))',squeeze(A(:,:,jr)));
end


