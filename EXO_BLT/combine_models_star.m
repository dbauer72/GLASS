function th_GL = combine_models_star(thc,thst);
% combine_models_star takes the marginal model for the star variables and
% the conditional model for the regional model and combines them into a
% larger model. 
%
% SYNTAX: th_GL = combine_models_star(thc,thst);
%
% INPUT: thc ... conditional model
%        thst ... marginal model
%
% OUTPUT: th_GL ... combined model.
%
% REMARK: star model comes first (otherwise we may have problems with
% echelon canonical form for stationary system). 
%
% AUTHOR: dbauer, 4.7.2025.

% calculate dimensions 
[ni,si] = size(thc.K);
[nist,sist] = size(thst.K);

% set up combined matrices 
Cf = zeros(si+sist,ni+nist);
Af = zeros(ni+nist,ni+nist);
Kf = zeros(ni+nist,si+sist);
Omegaf = zeros(si+sist,si+sist);


% calculate Dfull matrix to premultiply.
Df = eye(si+sist);
Df((sist+1):end,1:sist) = -thc.D;
iDf = inv(Df);

% calculate stacked C
Cf(1:sist,1:nist)=thst.C;
Cf((sist+1):end,(nist+1):end) = thc.C;

Cf = iDf*Cf;

% calculate stacked K
Kf(1:nist,1:sist)=thst.K;
Kf((nist+1):end,(sist+1):end) = thc.K;
Kf((nist+1):end,1:sist) = thc.B;

Kf = Kf*Df;
% calcuate stacked A
Af(1:nist,1:nist)=thst.A;
Af((nist+1):end,(nist+1):end) = thc.A;
Af((nist+1):end,1:nist) = thc.B*thst.C;

% calculate Omega
Omegaf(1:sist,1:sist)=thst.Omega;
Omegaf((sist+1):end,(sist+1):end) = thc.Omega;

Omegaf = iDf*Omegaf*iDf';

% fill matrices into theta object. 
th_GL = theta_urs;
th_GL.which = 'SS';
th_GL.A = Af;
th_GL.K=Kf;
th_GL.C=Cf;
th_GL.Omega = Omegaf; 

