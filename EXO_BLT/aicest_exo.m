function [k,Omega,AICs,exo_test,tharx,tharj,tharst] = aicest_exo(z,si,nmax,loo);
%
% 	estimates the lag length of a VARx model testing for exogeneity of the
% 	X-variables
%
% SYNTAX: [k,Omega,AICs,exo_test,tharx,tharj,tharst] = aicest_exo(z,sis,nmax);
%
% INPUTS:    z ... Tx (sist+si); [y* y], observations output (y) and input (y*)
%		     si ... dim(y)
%		  nmax ... maximal lag length
%          loo ... binary indicator, whether leave-one-out method should be
%                  used. 
%
% OUTPUTS:  k   ...  estimated lag order (joint system) according to AIC
%		    Omega ...  estimated innovation matrix. 
%           AICs ... vector of AIC values 
%         exo_test  ... structure: .Wald test statistic for exogeneity (Granger causality;
%                   leave one out method); .AIC values. 
%          tharx ... theta structure corresponding to ARX(k) 
%          tharj ... joint model
%          tharst ... AR(k) for star variables. 
%
% REMARKS: Within GLASS star variables y* are modelled as strongly
% exogenous (weak exogeneity plus Granger non-causality). This corresponds
% to the AR coefficients to be upper triangular. 
% This procedure estimates such a model and compares it to an unrestricted
% estimate with the same lag order. 
%
%  AUTHOR: dbauer, 14.7.2025

[T,nz] = size(z);
sist = nz - si;

if (nmax<0)
    k_est = -nmax;
    nmax = -nmax;
else
    k_est = [];
end

Aics = zeros(3,nmax+1);

Teff = T-nmax; % effective sample size.
Xfull = zeros(Teff,nmax*nz);
Xst = zeros(Teff,nmax*sist);
y = z(nmax+[1:Teff],sist+[1:si]);
yst = z(nmax+[1:Teff],1:sist);

for i=1:nmax
    Xfull(:,(i-1)*nz+[1:nz]) = z(nmax-i+[1:Teff],:);
    Xst(:,(i-1)*sist+[1:sist]) = z(nmax-i+[1:Teff],1:sist);
end


Xcond = [yst,Xfull];
% initialize: empty model. 
AICs(1,1) = log(det(yst'*yst/Teff))*Teff + Teff*(sist);
AICs(2,1) = log(det(y'*y/Teff))*Teff + Teff*(si);
AICs(3,1) = log(det(z'*z/Teff))*Teff + Teff*(sist+si);

for i=1:nmax
   
    % marginal model for star variables 
    res_ist = yst - Xst(:,1:sist*i)*(Xst(:,1:sist*i)\yst);
    Om_ist = res_ist'*res_ist/Teff;
    dOm_ist = det(Om_ist);
    if dOm_ist<10^(-20)
        dOm_ist = 10^(-20);
    end;
    AICs(1,i+1) = log(dOm_ist)*Teff + Teff*sist + 2*i*sist*sist;

    % conditional model 
    res_i = y - Xcond(:,1:nz*i)*(Xcond(:,1:nz*i)\y);
    Om_i = res_i'*res_i/Teff;
    dOm_i = det(Om_i);
    if dOm_i<10^(-20)
        dOm_i = 10^(-20);
    end;
    AICs(2,i+1) = log(dOm_i)*Teff + Teff*si  + 2*i*si*nz;

    % joint model 
    res_ij = [yst,y]  - Xfull(:,1:nz*i)*(Xfull(:,1:nz*i)\[yst,y]);
    Om_ij = res_ij'*res_ij/Teff;
    dOm_ij = det(Om_ij);
    if dOm_ij<10^(-20)
        dOm_ij = 10^(-20);
    end;
    AICs(3,i+1) = log(dOm_ij)*Teff + Teff*(sist+si)  + 2*i*nz*nz;
end;

% find optimal lag lengths
sig = AICs(1,:);
k(1) = find(~(sig > min(sig)))-1;

sig = AICs(2,:);
k(2) = find(~(sig > min(sig)))-1;

sig = AICs(3,:);
k(3) = find(~(sig > min(sig)))-1;

if isempty(k_est) % return the lag length selected for the conditional model. 
    k_est = k(2);
end;

%%%%% star model 
k_est = k(1);

phi = (Xst(:,1:(sist*k_est))\yst)';


if k_est>0
    if nargout>5
        tharst = theta_urs();
        tharst.which = 'poly'; %use polynomial form.
        tharst.a = [eye(sist),-phi];
        tharst.b = eye(sist);
        res = yst - Xst(:,1:(sist*k_est))*phi'; 
        tharst.Omega = res'*res/Teff;
        tharst.m = sist;
        tharst.num_param = length(phi(:));
    end
else
    tharst.a = eye(sist);
    tharst.b = eye(sist);
    tharst.Omega = yst'*yst/Teff;
    tharst.m = 0;
    tharst.num_param = 0;
end


%%%%% conditional model.
% fill in the model into the ARX theta structure for the conditional model 
k_est = k(2);
phi = (Xcond(:,1:(sist+nz*k_est))\y)';
res_i = y - Xcond(:,1:(sist+nz*k_est))*phi';
Omega = res_i'*res_i/Teff; 

if k_est>0
    if nargout>5
        tharx = theta_urs();
        tharx.which = 'poly'; %use polynomial form.
        % find indices for y and z terms
        indy = 2*sist+[1:si];
        if sist>0
            indz = [1:2*sist];
        else
            indz = [];
        end
        for j=2:k_est
            indy = [2*sist+[1:si],indy+nz];
            if sist>0
                indz = [indz,sist+(j-1)*nz+[1:sist]];
            end
        end;
        tharx.a = [eye(si),-phi(:,indy)];
        tharx.d = phi(:,indz);
        tharx.b = eye(si);
        
        tharx.Omega = Omega;
        tharx.m = sist;
        tharx.num_param = length(phi(:));
    end
else
    tharx.a = eye(si);
    tharx.d = zeros(si,sist);
    tharx.b = eye(si);
    tharx.Omega = y'*y/Teff;
    tharx.m = sist;
    tharx.num_param = 0;
end

%%%%% joint model 
k_est = k(3);

phi = (Xfull(:,1:(nz*k_est))\[yst,y])';
res_i = [yst,y]- Xfull(:,1:(nz*k_est))*phi';
Omega_j = res_i'*res_i/Teff;

if k_est>0
    if nargout>5
        tharj = theta_urs();
        tharj.which = 'poly'; %use polynomial form.
        % find indices for y and z terms
        indst = 1:(nz*k_est);
        tharj.a = [eye(nz),-phi(:,indst)];
        tharj.b = eye(nz);
        %res_i = [yst,y]- Xfull(:,1:(nz*k_est))*phi';
        tharj.Omega = Omega_j;
        tharj.m = 0;
        tharj.num_param = length(phi(:));
    end
else
    tharj.a = eye(sist);
    tharj.b = eye(sist);
    tharj.Omega = yst'*yst/Teff;
    tharj.m = 0;
    tharj.num_param = 0;
end

% exo test statistics
exo_test.AIC.vals = [AICs(3,k(3)+1),AICs(1,k(1)+1)+AICs(2,k(2)+1)]; 
exo_test.AIC.diff = AICs(3,k(3)+1)-AICs(1,k(1)+1)-AICs(2,k(2)+1); 
if (exo_test.AIC.diff>0)
    exo_test.AIC.dec = 'AIC: triangular model preferred.';
else
    exo_test.AIC.dec = 'AIC: unrestricted model preferred.';
end

% find the upper triangular part. 
if k(3)>0 % only if dynamics are contained. 
    phij= zeros(sist,(k_est-1)*si);
    indj = zeros(1,(k_est-1)*si);
    for i=1:(k_est-1)
        phij(:,(i-1)*si+[1:si]) = phi(1:sist,(i-1)*nz+sist+[1:si]);
        indj(1,(i-1)*si+[1:si]) = (i-1)*nz+sist+[1:si];
    end

    if (loo~=1) % surpluslag approach leaves out the last submatrix. 
        phij(:,(k_est-1)*si+[1:si]) = phi(1:sist,(k_est-1)*nz+sist+[1:si]);
        indj(1,(k_est-1)*si+[1:si]) = (k_est-1)*nz+sist+[1:si];
    end

    vPhi = phij(:);
    Omegar = tharj.Omega(1:sist,1:sist);
    iXtX = inv(Xfull(:,1:(nz*k_est))'*Xfull(:,1:(nz*k_est)));
    Vvphi = kron(inv(iXtX(indj,indj)),inv(Omegar));
    exo_test.Wald.val = vPhi'*Vvphi*vPhi;
    exo_test.Wald.df = length(vPhi);
    %exo_test.Wald.approx = (vPhi'*Vvphi*vPhi-exo_test.Wald.df)/sqrt(2*exo_test.Wald.df);
else
    exo_test.Wald.val = 0;
    exo_test.Wald.df = 0;
end; 





