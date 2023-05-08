function [V_theta,V_syst,T_syst_theta] = cal_quasi_like_V_X(param,y,s,m,n,c,D0,restrict);
% calculates the variance matrix for theta structure theta and the corresponding system.
%
% SYNTAX: [V_theta,V_syst,T_syst_theta] = cal_quasi_like_V_X(param,y,s,m,n,c,D0,Pbullet,restrict);
%
% INPUTS:  param ... d x 1 vector of parameter values (fed into param2syst,
%                       see there for description)
%          y     ... T x (s+m) matrix of observations of endo- and
%                    exogenous variables;
%          s     ... integer; dimension of endogenous vars
%          m     ... integer; dimension of exogenous vars (e.g.
%                     deterministics)
%          n     ... integer; state dimension
%          c     ... integer; number of common trends.
%          D0  ... indicator; If D0=1 (TRUE), then D=0 is used.
%          restrict ... structure containing restrictions to be passed on
%                       to param2syst. 
%
% OUTPUT:   V_theta ... real matrix; variance of parameter estimators
%           V_syst  ... real matrix; variance of system matrix estimators. 
%           T_syst_theta ... real matrix containing the derivatives of the
%           vectorizations of the system matrices as a function of the
%           parameters. 
%
% REMARKS: 
%         + the first parameters corresponding to Omega are disregarded. 
%         + restrict.Omega: if this field is included, Omega is fixed not
%                         included in parameter vector.
%         + not implemented yet!
%         
% AUTHOR: dbauer, 27.4.2023

if nargin<8
    restrict.det_res =0;
    restrict.scale = ones(length(param),1);
end;

if isfield(restrict,'Omega')
    Om = restrict.Omega;
    paraomi = extr_lowtri(Om);
    param = [paraomi(:);param];
end;

% extract matrices from parameters 
[A,B,C,D,K,Omega,~,gr_syst] = param2syst_X(param,s,n,m,c,D0,restrict); 

nparOm = length(extr_lowtri(Omega));

% convert params to matrices
u = y(:,s+1:end);
y = y(:,1:s);

T = size(y,1);

Abar = A-K*C;
xh = ltitr(Abar,[B-K*D,K],[u,y],zeros(n,1));
tres = y-xh*C'-u*D';
[tres,x_0] = est_initial_val(tres,Abar,K,C);
xh = ltitr(Abar,[B-K*D,K],[u,y],x_0);

% initialize matrices 
npar = length(gr_syst)-nparOm;
gr_syst = gr_syst(nparOm+1:end); 

V_theta = zeros(npar,npar); 

nentries = n^2 + 2*s*n + s*m; % A, K, C, D. 
T_syst_theta = zeros(nentries,npar); 
V_syst = zeros(nentries,nentries);

% cycle over entries of parameters
dres = zeros(size(tres,1),size(tres,2),npar); 

for a=1:npar
    % derivative trafo params to syst 
    T_syst_theta(1:n^2,a) = gr_syst(a).A(:);
    T_syst_theta(n^2+[1:n*s],a) = gr_syst(a).C(:);
    T_syst_theta(n^2+n*s+[1:n*s],a) = gr_syst(a).K(:);
    T_syst_theta(n^2+2*n*s+[1:m*s],a) = gr_syst(a).D(:);

    % derivative of eps 
    dAbar = gr_syst(a).A -gr_syst(a).K*gr_syst(a).C;
    dK = gr_syst(a).K;
    dC = gr_syst(a).C;
    dD = gr_syst(a).D;

    % input to state equation
    dKe = xh*dAbar' + y*dK' - u*dD'*K';

    dx = ltitr(Abar,eye(n),dKe,zeros(n,1));
    dres(:,:,a) = - u*dD' - xh*dC' - dx*C';
end;

iOM = inv(Omega);
H_theta = V_theta;
% now start the calculation of the variance
for a= 1:npar
    depsa = squeeze(dres(:,:,a));
    H_theta(a,a) = trace(iOM* depsa'*depsa); 
    for b=(a+1):npar
        % derivative of eps
        depsb = squeeze(dres(:,:,b));
        H_theta(a,b) = trace(iOM* depsa'*depsb); 
        H_theta(b,a) = H_theta(a,b);
    end
end
meH = min(eig(H_theta));

if (meH<0.00001) % eigenvalue is too small? 
    H_theta = H_theta + eye(npar)*(meH+0.0001);
end

V_theta = inv(H_theta);
V_syst = T_syst_theta * V_theta * T_syst_theta';
