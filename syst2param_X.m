function par = syst2param_X(th,c,restrict);
% syst2param_X collects the parameters corresponding to the system subject to
% the restrictions formulated in restrict. It also uses the matrices
% corresponding to exogenous inputs. 
%
% SYNTAX: par = syst2param_X(th,c,restrict);
%
% INPUT: th ... theta structure
%        c  ... number of common trends
%        restrict ... hypotheses 
% 
% OUTPUT: par ... vector of parameters.
%
% REMARK: currently restricted to stationary case!! c and restrict are
% inactive. 
%
% AUTHOR: dbauer, 27.4.2023. 

paraomi = extr_lowtri(th.Omega);
parami = th2param_X(th,c,1);

par = [paraomi(:);parami(:)];

%if (restrict.det_res)&&(c>0) % restrict the highest order deterministic term.    
%    C1 = th.C(:,1:c);
%    D = th.D; 
%    pard = C1'*D(:,1);
%    D2 = D(:,2:end);
%    if size(D,2)>1
%        pard = [pard(:);D2(:)];
%    end;
%    par = [paraomi(:);parami(:);pard];
%else
%    par = [paraomi(:);parami(:);th.D(:)];
%end