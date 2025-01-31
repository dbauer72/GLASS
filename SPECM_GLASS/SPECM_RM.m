function  [result,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM_RM(z,n,index,nmax,Pbull,thi);
% SPECM_RM optimizes the Gaussian pseudo likelihood starting from an inital
% CCA estimated system for the joint regional variables in the context of GLASS.
% The star and the dominant variables are included in the modelling
% resulting in a larger dimensional system for the i-th regional model. 
% 
% SYNTAX: [result,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM_RM(z,n,index,nmax,Pbull,thi);
%
% INPUT:  z ... Tx(si+sist) data: [y_{i,t},y_{i,t}^*]
%         n ... state order to use
%         index ... vector of integers; see 'cal_quasi_like_RM' for a
%                   description
%         nmax ... maximum for AR order selection.
%         Pbull ... indicator; if Pbull>0, state is started in stationary
%                   distribution.
%         thi    ... use initial estimate rather than CCA.
%
% OUTPUT: result ... structure containing the estimation results.
%         th     ... theta estimate
%         (Ac,Kc,Cc) ... state space system in canonical form.
%         Omegac ... noise variance estimate
%         thi   ... initial estimate corresponding to CCA.
%         thi2  ... initial estimate after RH-procedure.
%         lle   ... minimizing value of the unrestricted quasi scaled log
%                  likelihood (over M_n). 
%
% AUTHOR: dbauer, 2.8.2024.

% define restrict to comply with functions from EICIP. 
restrict.det_res = 0;

if nargin<5
    Pbull = -1;
end;

[T,sf]= size(z); 
si = index(1);
ci = index(2);
cist = index(3);

if ci>n
    disp('Wrong specification: c cannot be larger than n!');
    result = est_result();
    return;
end;
if ci<cist
    disp('Wrong specification: c^* cannot be larger than c!');
    result = est_result();
    return;
end;
if sf<si
    disp('Wrong specification: sf cannot be larger than s!');
    result = est_result();
    return;
end;


% STEP 1: estimate an unrestricted model using the procedures for
% stationary data with no exogenous inputs/deterministics present. 
%
if nargin<6
    % initial VAR estimate
    [k,sig,kbc,sigbc,phi]  = aicest(z,size(z,2),nmax);
    
    % CCA estimate
    if 2*k>sqrt(T) % reduce k if too big. 
        k = max(1,floor(sqrt(T)/2)); 
    end
    if 2*k*sf<n
        k = n;
    end;
    [thi,Ai,Ki,Ci,Omegai] = CCA(z,n,2*k,2*k,0);
end;

% find parameters in stationary version
parami = syst2param(thi,0,restrict);

% improve estimate
options = optimoptions('fminunc','display','final');
options.MaxFunctionEvaluations = 10000;
restrict.scale = ones(length(parami),1);

pare = est_cal_like_hess(z,sf,0,n,0,thi,restrict);

% random improvement of stable systems
result = compile_results(pare,sf,0,n,0,restrict,z(:,1:sf),zeros(T,0),Pbull,nmax);
lle_c = result.deviance;
lle=lle_c+1;
it = 0;
while ((it<5)&&(lle_c<lle))
    if it>0
        disp(sprintf('Running random improve %2d. Old value: %4.2f -> new value: %4.2f.',it,lle,lle_c));
    end
    lle = lle_c;
    result = random_improve_call(result,0.001,5);
    lle_c =result.deviance;
    it = it+1;
end;
pare = result.param;
lle =result.deviance; %cal_quasi_like(pare,z,s,m,n,0,Pbull,restrict);

[~,~,~,~,~,the] = param2syst(pare,sf,n,0,0,restrict);
thi2=the; % save the intermediate result.

% STEP 2: based on the initial estimate find closest state space system with c common trends, if c>0. 
% Else skip this step. 
if ci>0 
    [LL,alphahat,betahat,Chat] = RH_specm(z(:,1:sf),the.A-the.K*the.C,the.K,sf-ci,'y');

    Abar = the.A-the.K*the.C;
    thi2=the;
    the.A =Abar+the.K*Chat;
    the.C = Chat;

    % find initial estimate.
    % find corresponding params
    parami = syst2param(the,ci,restrict);
    restrict.scale = ones(length(parami),1);

    parc = est_cal_like_hess(z,sf,0,n,ci,the,restrict);
    [~,~,~,~,~,the] = param2syst(parc,sf,n,0,ci,restrict);

% STEP 2b: change the parameterization and approximate the K_{i,c}^* matrix,
% afterwards optimize again. 
% 
% Here comes the GLASS specific part!!
    parc_star = th2param_RM(the,index,1);
    paromi = extr_lowtri(the.Omega);
    param_star = [paromi(:);parc_star(:)];
    [pare_star,fval,exitflag] = fminunc(@(x) cal_quasi_like_RM(x,z,n,index,Pbull),param_star,options);
    parc = pare_star;
else % if c==0 (stationary case)
    parc = pare;
end

% STEP 4: collect results 
result = compile_results_RM(parc,n,index,z(:,1:sf),Pbull,nmax);

% STEP 5: convert into regional vars and exogenous vars 
thc = result.theta;
Ac = thc.A;
Kc = thc.K;
Cc = thc.C;
Omegac = thc.Omega;
