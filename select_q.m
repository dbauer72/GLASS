function q = select_q(Omega,T); 
% select_q returns the value of the integer q estimated as the rank of the
% matrix Omega using information criteria.
%
% SYNTAX: q = select_q(Omega);
% 
% INPUT:  Omega rxr real covariance matrix. 
%
% OUTPUT: q ... integer.
%
% AUTHOR: dbauer, 30.11.2023. 

r = size(Omega,1);

[U,S,V]= svd(Omega);
 
% q obtained by thresholding. Everything smaller than Q_T is discarded. 
QT= log(T)/sqrt(T);
q = sum(diag(S) >QT); 
