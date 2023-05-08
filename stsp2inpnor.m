function th = stsp2inpnor(th);
% converts the singular state space system to an input normal state space system.
% 
%  SYNTAX: th = stsp2inpnor(th);
%  
% INPUT:  th .. theta structure.
%
% OUTPUT: th .. theta structure.
%
% REMARK: uses the system (A,B,C,D) and converts it to a form such that 
% + P = eye(n): State variance is the identity matrix.
% + Q = diag(S): observation Gramian is diagonal with sorted SVs. 
% + Columns of C start with a positive entry. 
% + norm(D) = 1 

A = th.A;
B = th.B;
C = th.C;
D = th.D; 
Om = th.Omega;
pD = pinv(D); 

[N,n]= size(C);
q = size(D,2); 

% first step: input normal. 
BOB = B*th.Omega*B';
P = reshape(inv(eye(n^2) - kron(A,A))*BOB(:),n,n);
[U,S] = svd(P); 
Tr = U*diag(diag(sqrt(S)));
iTr = inv(Tr); 

A = iTr*A*Tr;
B = iTr*B;
C = C*Tr;

% second step: diagonalize observability Gramian. 
CtC = C'*C;

Q = reshape(inv(eye(n^2) - kron(A',A'))*CtC(:),n,n);
[U,S]= svd(Q);
Tr = U;
iTr = U';
A = iTr*A*Tr;
B = iTr*B;
C = C*Tr;

% third step: first entries in C positive. 
si = sign(C(1,:));
Tr = diag(si);
iTr = Tr; 
A = iTr*A*Tr;
B = iTr*B;
C = C*Tr;

% last step: normalize k(1)*Omega*k(1)'.
kv1 = D + C*inv(eye(n)-A)*B;
DOD = kv1*Om*kv1';
[U,S]= svd(DOD); 

kv1n = U(:,1:q);
Om = S(1:q,1:q); 

si = sign(kv1n(1,:));
kv1n = kv1n*diag(si);
Tr = kv1\kv1n;
B = B*Tr;
D = D*Tr; 

% fill in normalized matrices.
th.A = A;
th.B = B;
th.C = C;
th.D = D;
th.Omega = Om; 

