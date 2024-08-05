function result = compile_results_RM(parc,n,index,y,Pbull,nmax);
% compile_results_RM writes all information corresponding to the parameter
% vector parc in the i-th regional model into the structure result. 
% 
% SYNTAX: result = compile_results(parc,n,index,y,Pbull,nmax);
% 
% INPUT: parc ... parameter vector; is passed tp param2syst.
%        n    ... integer; state dimension
%        index    ... vector of integers; see 'cal_quasi_like_RM' for
%                   detailed explanation
%        y    ... Txs; observations of endogenous vars.
%        Pbull ... indicator; if Pbull>0: state started with stationary
%                  dist; else with zero.
%       nmax  ... maximal state dimension allowed (contained for
%       information in result.call. 
%
% OUTPUT: result: est_result structure (see there for description of
% fields.
%
% REMARKS: + system is calculated using param2th_RM; 
%          + quasi likelihood calculated using cal_quasi_like_RM. 
%       
% AUTHOR: dbauer, 2.8.2024.

[T,sf]=size(y);
[thc,paromi] = param2th_RM(parc,n,sf,index); 

if length(paromi)>0
    Omega = fill_lowtri(paromi,sf);
else
    Omega = eye(sf); 
end
thc.B = zeros(n,m);

[llc,resc] = cal_quasi_like_RM(parc,y,n,index,Pbull);
Omega= resc'*resc/size(resc,1);
thc.Omega = Omega;

% calculate the variance matrices 
%[V_theta,V_syst,~] = cal_quasi_like_V(parc,[y,dt],s,m,n,c,restrict);

% fill in everything into result structure
result = est_result();
result.urs = index;
if index(2)>0
    result.ur = 'I(1)';
else
    result.ur = 'I(0)';
end

result.n = n;
result.s = index(1);
result.y = y;
result.dt = [];
result.theta = thc;
result.deviance = llc;
result.aic = llc + 2*length(parc);
result.bic = llc + log(size(resc,1))*length(parc);
result.res = resc;
result.restrict = [];
result.Pbull = Pbull;
result.call = sprintf('SPECM_RM(z,%d,[%d,%d,%d],%d,%d)',n,index,nmax,Pbull);
result.param = parc;
%result.V_theta = V_theta;
%result.V_syst = V_syst;
