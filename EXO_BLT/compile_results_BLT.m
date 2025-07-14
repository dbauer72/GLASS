function result = compile_results_BLT(parc,n,sf,index,y,Pbull,nmax,restrict);
% compile_results_RM writes all information corresponding to the parameter
% vector parc in the i-th regional model into the structure result. 
% 
% SYNTAX: result = compile_results_BLT(parc,n,sf,index,y,Pbull,nmax);
% 
% INPUT: parc ... parameter vector; is passed tp param2syst.
%        n    ... integer; state dimension
%        sf   ... pair of integers; dimension vars region i; dim ystar. 
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
% AUTHOR: dbauer, 4.6.2025
[T,stot]=size(y);
si = sf(1);
sist = sf(2);
ci = index(1);
cist = index(2); 

% extract matrices from parameters
if nargin<8
    restrict.det_res = 0;
end

[thc] = param2th_BLT(parc,n,sf,index); 

[llc,resc] = cal_quasi_like_BLT(parc,y,n,sf,index,Pbull,restrict);
Omega= resc'*resc/size(resc,1);
thc.Omega = Omega;

% calculate the variance matrices 
%[V_theta,V_syst,~] = cal_quasi_like_V(parc,[y,dt],s,m,n,c,restrict);

% fill in everything into result structure
result = est_result();
result.urs = index;
if index(1)>0
    result.ur = 'I(1)';
else
    result.ur = 'I(0)';
end

result.n = n;
result.s = sf;
result.y = y;
result.dt = [];
result.theta = thc;
result.deviance = llc;
result.aic = llc + 2*length(parc);
result.bic = llc + log(size(resc,1))*length(parc);
result.res = resc;
result.restrict = index;
result.Pbull = Pbull;
result.call = sprintf('SPECM_BLT(z,%d,[%d,%d,%d,%d],%d,%d)',n,sf,index,nmax,Pbull);
result.param = parc;
%result.V_theta = V_theta;
%result.V_syst = V_syst;
