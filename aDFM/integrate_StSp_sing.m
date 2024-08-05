function ith = integrate_StSp_sing(th);
% provides the transfer function corresponding to the integrated series. 
%
% SYNTAX: dth = integrate_StSp_sing(th);
%
% INPUT: th ... theta structure containing the transfer function in state
% space form assumed to be singular. 
%  
% OUTPUT: ith ... theta structure for integrated process. 
[n,q]= size(th.B);
r = size(th.C,1);

ith =th;

A= [th.A,zeros(n,r);th.C,eye(r)];
B = [th.B;th.D];
C = [th.C,eye(r)];

% transformation 
Trafo = [eye(n),zeros(n,r);-th.C*inv(th.A-eye(n)),eye(r)];

tA = Trafo*A*inv(Trafo);
tB = Trafo*B;
tC = C*inv(Trafo);

% reduce order
B2 = tB(n+[1:r],:);
[u,s,v] = svd(B2); 

% fill in reduced transfer function 
ith.A = [tA(1:n,1:n),zeros(n,q);zeros(q,n),eye(q)];
ith.B = [tB(1:n,:);s(1:q,1:q)*v(:,1:q)'];
ith.C = [tC(:,1:n),tC(:,n+[1:r])*u(:,1:q)];



