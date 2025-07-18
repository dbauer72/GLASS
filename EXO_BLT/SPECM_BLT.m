function  [result,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM_BLT(z,n,sf,index,nmax,Pbull,thi);
% SPECM_BLT optimizes the Gaussian pseudo likelihood starting from an inital
% CCA_X estimated system for the joint regional variables in the context of GLASS.
% The star and the dominant variables are included in the modelling
% resulting in a larger dimensional system for the i-th regional model. 
% 
% SYNTAX: [result,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM_BLT(z,n,sf,index,nmax,Pbull,thi);
%
% INPUT:  z ... Tx(si+sist) data: [y_{i,t},y_{i,t}^*]
%         n ... state order to use
%        sf   ... pair of integers; dimension vars region i; dim ystar. 
%         index ... vector of integers; see 'cal_quasi_like_RM' for a
%                   description
%         nmax ... maximum for AR order selection.
%         Pbull ... indicator; if Pbull>0, state is started in stationary
%                   distribution.
%         thi    ... use initial estimate rather than CCA_X.
%
% OUTPUT: result ... structure containing the estimation results.
%         th     ... theta estimate
%         (Ac,Bc,Cc,Dc,Kc) ... state space system in canonical form.
%         Omegac ... noise variance estimate
%         thi   ... initial estimate corresponding to CCA_X.
%         thi2  ... initial estimate after RH-procedure (for I(1) syst).
%         lle   ... minimizing value of the unrestricted quasi scaled log
%                  likelihood (over M_n). 
%
% AUTHOR: dbauer, 4.6.2025.

% define restrict to comply with functions from EICIP. 
restrict.det_res = 0;

if nargin<6
    Pbull = 0;
end;

[T,stot]= size(z); 
si = sf(1);
sist = sf(2);
ci = index(1);
cist = index(2); 


if ci>n
    disp('Wrong specification: c cannot be larger than n!');
    result = est_result();
    return;
end;
if (stot ~= sum(sf))
    disp('Wrong specification: number of columns in y must equal sum(sf)!');
    result = est_result();
    return;
end;


% STEP 1: estimate an unrestricted model using the procedures for
% stationary data with no exogenous inputs/deterministics present. 
%
if nargin<6
    % initial VAR estimate
    [k,sig,kbc,sigbc,phi]  = aicest(z,si,nmax);
    
    % CCA estimate
    if 2*k>sqrt(T) % reduce k if too big. 
        k = max(1,floor(sqrt(T)/2)); 
    end
    if 2*k*sf<n
        k = n;
    end;
    [thi,Ai,Bi,Ci,Di,Ki,Omegai] = CCAX(z(:,[(sist+1):end,1:sist]),si,n,2*k,2*k,0,1);
    parami = th2param_BLT(thi,[0,index(2)],1);
    paro = extr_lowtri(Omegai);
    [qlikei,tresi] = cal_quasi_like_BLT([parami(:);paro(:)],z,n,[si,sist],[0,index(2)],0);

    % alternative: long VARX and project
    [k,Omega,AICs,exo_test,tharx,thj,thst] = aicest_exo(z,si,2*k,1);
    thi_AR = project_init_ARX(tharx,si,sist,n,n,0,0);
    parami2 = th2param_BLT(thi_AR,[0,index(2)],1);
    paro2 = extr_lowtri(thi_AR.Omega);
    [qlikei2,tresi2] = cal_quasi_like_BLT([parami2(:);paro2(:)],z,n,[si,sist],[0,index(2)],Pbull);

    if (qlikei2<qlikei)
        thi=thi_AR;
    end

end;

% find parameters in stationary version
parami = th2param_BLT(thi,[0,index(2)],1);

% improve estimate
options = optimoptions('fminunc','display','final');
options.MaxFunctionEvaluations = 10000;
%restrict.scale = ones(length(parami),1);
restrict.det_res = 0; 

% estimate initial stable model
pare = est_cal_like_hess_BLT(z,n,sf,[0,index(2)],thi,Pbull,restrict);

% random improvement of stable systems
result = compile_results_BLT(pare,n,sf,[0,index(2)],z,Pbull,nmax);
lle_c = result.deviance;
lle=lle_c+1;
it = 0;
while ((it<5)&&(lle_c<lle))
    if it>0
        disp(sprintf('Running random improve %2d. Old value: %4.2f -> new value: %4.2f.',it,lle,lle_c));
    end
    lle = lle_c;
    result = random_improve_BLT(result,0.00001,5);
    lle_c =result.deviance;
    it = it+1;
end;
pare = result.param;
lle =result.deviance; %cal_quasi_like(pare,z,s,m,n,0,Pbull,restrict);

[the] = param2th_BLT(pare,n,sf,[0,index(2)]);
thi2=the; % save the intermediate result.

% STEP 2: based on the initial estimate find closest state space system with c common trends, if c>0. 
% Else skip this step. 
if ci>0 
    [LL,alphahat,betahat,Chat] = RH_specm_BLT(z,the.A-the.K*the.C,the.B-the.K*the.D,the.K,cist,si-ci,'y');

    Abar = the.A-the.K*the.C;
    thi2=the;
    the.A =Abar+the.K*Chat;
    the.C = Chat;

    % find initial estimate.
    % find corresponding params
    parami = th2param_BLT(the,index);
    %restrict.scale = ones(length(parami),1);

    parc = est_cal_like_hess_BLT(z,n,sf,index,the,Pbull,restrict);

else % if c==0 (stationary case)
    parc = pare;
end

% STEP 4: collect results 
result = compile_results_BLT(parc,n,sf,index,z,Pbull,nmax);

thc = result.theta;
Ac = thc.A;
Bc = thc.B;
Dc = thc.D;
Kc = thc.K;
Cc = thc.C;
Omegac = thc.Omega;
