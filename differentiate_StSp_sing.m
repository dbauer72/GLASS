function dth = differentiate_StSp_sing(th);
% provides the transfer function corresponding to the differenced series. 
%
% SYNTAX: dth = differentiate_StSp_sing(th);
%
% INPUT: th ... theta structure containing the transfer function in state
% space form assumed to be singular. 
%  
% OUTPUT: dth ... theta structure for differenced process. 
[n,q]= size(th.B);
r = size(th.C,1);

dth =th;

dth.A= [th.A,th.B;zeros(q,n+q)];
dth.B = [zeros(n,q);eye(q)];
dth.C = [th.C*(th.A-eye(n)),th.C*th.B-th.D];
