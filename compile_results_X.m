function result = compile_results_X(parc,s,m,n,c,restrict,y,u,D0,Pbull,nmax);
% compile_results writes all information corresponding to the parameter
% vector parc into the structure result. 
% 
% SYNTAX: result = compile_results(parc,s,m,n,c,restrict,y,dt,Pbull,nmax);
% 
% INPUT: parc ... parameter vector; is passed tp param2syst.
%        s    ... integer; dimension of endogenous vars
%        m    ... integer; dimension of exogenous vars
%        n    ... integer; state dimension
%        c    ... integer; number of common trends
%        restrict ... structure; contains information for parameterization.
%        y    ... Txs; observations of endogenous vars.
%        u   ... Txm; observations of exogenous vars.
%        D0  ... indicator; If D0=1 (TRUE), then D=0 is used. 
%        Pbull ... indicator; if Pbull>0: state started with stationary
%                  dist; else with zero.
%       nmax  ... maximal state dimension allowed (contained for
%       information in result.call. 
%
% OUTPUT: result: est_result structure (see there for description of
% fields.
%
% REMARKS: + system is calculated using param2syst_X; 
%          + quasi likelihood calculated using cal_quasi_like_X.
%          + implements only stationary case for now. 
%       
% AUTHOR: dbauer, 27.4.2023.
if isfield(restrict,'scale')
    restrict = rmfield(restrict,'scale');
end;
[~,~,~,~,~,~,thc] = param2syst_X(parc,s,n,m,c,D0,restrict);


if Pbull<0 
    sizOm = s*(s+1)/2;
else
    sizOm = 0;
end;
[llc,resc] = cal_quasi_like_X(parc(sizOm+1:end),[y,u],s,m,n,c,D0,Pbull,restrict);
Omega= resc'*resc/size(resc,1);
th.Omega = Omega;

% calculate the variance matrices 
[V_theta,V_syst,~] = cal_quasi_like_V_X(parc,[y,u],s,m,n,c,D0,restrict);

% fill in everything into result structure
result = est_result();
result.urs = c;
if c>0
    result.ur = 'I(1)';
else
    result.ur = 'I(0)';
end

result.n = n;
result.s = s;
result.y = y;
result.dt = u;
result.theta = thc;
result.deviance = llc;
result.aic = llc + 2*length(parc);
result.bic = llc + log(size(resc,1))*length(parc);
result.res = resc;
result.restrict = restrict;
result.D0 = D0;
result.Pbull = Pbull;
result.call = sprintf('SPECM(z,%d,%d,%d,%d,%d)',s,n,c,nmax,Pbull);
result.param = parc;
result.V_theta = V_theta;
result.V_syst = V_syst;
