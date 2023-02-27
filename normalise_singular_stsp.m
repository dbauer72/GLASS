function the = normalise_singular_stsp(the);
% normalise_singular_stsp provides a system transformation such that the
% noise variance is the identity, C is orthnormal, Omega = I, D is p.l.t. 
% and the controllability matrix is p.u.t. 
%
% SYNTAX:  the = normalise_singular_stsp(the);
%
% INPUT: the 
%
% OUTPUT: the
%
% AUTHOR: dbauer 27.1.2023. 

% correct for Omega. 
q =size(the.Omega,2);
r = size(the.A,1); % state dimension
Sig = the.D* the.Omega * the.D';
[u,s,v] = svd(Sig);
De = u(:,1:q)*diag(sqrt(diag(s(1:q,1:q))));

[Q,R]=qr(De');
De = R'; 
Tr = the.D\De;
the.D = R';
the.Omega = eye(q); 
the.B = the.B*Tr; 

% make C orthonormal.
C = the.C;
[u,s,v]= svd(C);
the.C = u(:,1:r);
Tr = s(1:r,1:r)*v'; 
A = Tr*the.A*inv(Tr);
B = Tr*the.B(:,1:q);
% controllability 
Contr = zeros(r,r*q);
AjB = B;
for j=1:r
    Contr(:,(j-1)*q+[1:q]) = AjB;
    AjB = A*AjB;
end; 

[Q,R]= qr(Contr);
the.C = the.C*Q;
the.A = Q'*A*Q;
the.B = R(:,1:q);
