function [qlike,tres] = cal_quasi_like_RM(param,y,n,index,Pbullet,restrict);
% calculates the quasi likelihood for theta structure theta in the regional model within the GLASS framework.
%
% The structure is characterized by the index = [s(i),c(i),c(i)^*]. 
% Here s(i) denotes the dimension of the regional variables, 
%      c(i) the number of common trends in the regional model i
%      c(i)^* the rank of K_{i,c}^* indicating the common trends spilling
%      over. 
%
% SYNTAX: [qlike,tres] = cal_quasi_like_RM(param,y,s,m,n,c,Pbullet,restrict);
%
% INPUTS:  param ... d x 1 vector of parameter values (fed into param2syst,
%                       see there for description)
%          y     ... T x s matrix of observations
%          n     ... integer; state dimension
%          index ... vector of integers; see above
%          Pbullet ... indicator var; if >0: state started with stationary
%                   variance; if 0: state initialized with zero (corresponds to prediction error).  
%          restrict ... restrict.exo (if it is present implements
%                       exogeneity of star and dominant variables).
%
% OUTPUT:   qlike ... real; -2/T log Gaussian likelihood.
%           tres  ... Txs; matrix of residuals.
%           xts   ... Txn; matrix of state estimates.
%
% REMARKS: 
%         + large penalization of eigenvalues larger than 0.99 of tilde A (unit
%         roots) and tilde A-KC (minimum phase). 
%         
% AUTHOR: dbauer, 2.8.2024

if nargin<6
    restrict.det_res = 0;
end

if nargin<5 % Pbullet: start with x_1 = x_bullet? 
    Pbullet = 0;
end;

qlike = 0;
[T,sf] = size(y);
si = index(1);
ci = index(2);
cist = index(3); 

% extract matrices from parameters
if isfield(restrict,'exo')
    Omega= zeros(sf,sf);
    sist = sf-si; 
    sizOmi = si*(si+1)/2;
    paromi = param(1:sizOmi);
    Omega(1:si,1:si)= fill_lowtri(paromi,si);

    sizOmist = sist*(sist+1)/2;
    paromist = param(sizOmi+[1:sizOmist]);
    Omega(si+[1:sist],si+[1:sist])= fill_lowtri(paromist,sist);

    sizOm = sizOmi+sizOmist; 
else
    sizOm = sf*(sf+1)/2;
    paromi = param(1:sizOm);
    Omega = fill_lowtri(paromi,sf);
end
[th] = param2th_RM(param(sizOm+1:end),n,sf,index); 
A= th.A;
K= th.K;
C= th.C;

%if length(paromi)
%    Omega = fill_lowtri(paromi,sf);
%else
%    Omega = eye(sf); 
%end

% new switch for Pbullet == 0: PE estimation 
% --- if not minimum-phase-> penalize and project! ---
makcA = max(abs(eig(A-K*C)));
if (makcA>0.99) && (Pbullet<1)
    A = A*0.99/makcA;
    K = K*0.99/makcA;
    qlike = 10^6*(exp(makcA)-exp(0.99));
end


if Pbullet == 0 
    xh = ltitr(A-K*C,K,y,zeros(n,1));
    tres = y-xh*C';
    [tres,~] = est_initial_val(tres,A-K*C,K,C);
    Omegah = tres'*tres;
    qlike = qlike+ T*log(det(Omega)) + trace(inv(Omega)*Omegah);
    return;
end;

if Pbullet<0 
    xh = ltitr(A-K*C,K,y,zeros(n,1));
    tres = y-xh*C';
    [tres,~] = est_initial_val(tres,A-K*C,K,C);
%     x1e = zeros(T*s,n);
%     x1e([1:s],:)=C;
%     for j=2:T
%         x1e((j-1)*s+[1:s],:)= x1e((j-2)*s+[1:s],:)*(A-K*C);  
%     end
%     ttres = tres';
%     ttv = ttres(:);
%     ttv2 = ttv - x1e*(x1e\ttv);
%     tres = reshape(ttv2,s,T)';
    Omegah = tres'*tres/T;
    qlike = qlike+ T*log(det(Omegah)) + T*sf; %trace(inv(Omega)*Omegah);
    return;
end


if ci>0 % there is an integration going on! 
    % find tilde k. via C1: 
    C1 = C(:,1:ci);
    Pi = C1*inv(C1'*C1)*C1';
    Cbull = [C(:,ci+1:end)];
    Kbull= K(ci+1:end,:);
    K1 = K(1:ci,:);

    % calculate tilde y_t = y_t - Pi y_{t-1}
    ty = y; 
    ty(2:end,:)=ty(2:end,:)-y(1:end-1,:)*Pi;

    % transformed matrices 
    tilA = [zeros(ci,ci),-C1'*Cbull;A(ci+1:end,:)];
    tilK = [K1-C1';Kbull];
    tilC = C;
else % no integration
    tilA = A;
    tilK= K;
    tilC=C;
    ty = y;
end;

if Pbullet>0 % if P_bull to be included: system cannot be unstable!
    maA = max(abs(eig(tilA)));
    if maA >0.99
        tilA = tilA*0.99/maA;
        qlike = qlike + 10^6*(exp(maA)-exp(0.99));
    end
end
% correct for not-invertible systems: add large penalty.
makcA = max(abs(eig(tilA-tilK*tilC)));
if makcA>0.99
    tilA = tilA*0.99/makcA;
    tilK = tilK*0.99/makcA;
    qlike = qlike + 10^6*(exp(makcA)-exp(0.99));
end

% calculate P_bull
Q = tilK*Omega*tilK';
Q = (Q+Q')/2;
Qbull = Q(c+1:end,c+1:end);
if Pbullet>0
    pb = inv( eye((n-ci)^2) - kron(tilA(ci+1:end,ci+1:end),tilA(ci+1:end,ci+1:end)))*Qbull(:);
    Pbull = reshape(pb,n-ci,n-ci);
    P0 = [0*eye(ci),zeros(ci,n-ci);zeros(n-ci,ci),Pbull];
    
    
    P0 = (P0+P0')/2;
    % --- initialize the Kalman filter ---
    x0= zeros(n,1);
    tres = ty*0;
    tres(1,:)=ty(1,:)-x0'*C';
    Omegat= tilC*P0*tilC'+Omega;
    Kt = (tilA*P0*tilC'+tilK*Omega)*inv(Omegat);
    xf = tilA*x0 + Kt*tres(1,:)';
    Pkg1 = tilA*P0*tilA' + Q- Kt*Omegat*Kt'; % P(2|1)
    
    %xf = x0 + P0*tilC'*inv(tilC*P0*tilC'+Omega)*(y(1,:)'-tilC*x0);
    %P0g0 = P0-P0*tilC'*inv(tilC*P0*tilC'+Omega)*tilC*P0;
    %
    %% --- run the Kalman filter ---
    %
    %xf = tilA*xf;
    %Pkg1 = tilA*P0g0*tilA' + Q;
    %tres = ty*0;
    %tres(1,:)=ty(1,:);
    %Omegat= tilC*P0*tilC'+Omega;
    
    qlike = qlike +  log(det(Omegat)) + tres(1,:)*inv(Omegat)*tres(1,:)';
    
    for t=2:T % filter equations
        Omegat = tilC*Pkg1*tilC'+Omega;
        iOm = inv(Omegat);
        tres(t,:)= ty(t,:) - xf'*tilC';
        Kt = (tilA*Pkg1*tilC'+tilK*Omega)*iOm;
        xf = tilA*xf+ Kt*tres(t,:)';
        Pkg1 = (tilA * Pkg1 *tilA') +Q - (Kt*Omegat*Kt');
        Pkg1 = (Pkg1+Pkg1')/2;
        % update likelihood
        qlike = qlike + log(det(Omegat)) + tres(t,:)*inv(Omegat)*tres(t,:)';
    end
    
else
    P0 = 0*eye(n);
    x0= zeros(n,1); % x(1|0)
    Abar = tilA-tilK*tilC;
    % --- run the filter for the inverse  transfer function ---
    Ts = size(ty,1);
    xt = ltitr(Abar,tilK,ty(1:Ts-1,:),x0);
    tres = ty(1:Ts,:)-xt*tilC';
    [tres,~] = est_initial_val(tres,Abar,tilK,tilC);
    
    Omegat =  tres'*tres/(Ts);
    % update likelihood
    qlike = Ts*log(det(Omegat)) + sf*Ts;
end


%qlike = qlike/T;
