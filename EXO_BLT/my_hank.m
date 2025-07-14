function H = my_hank(A,K,C,krow,kcol)
% calculates the Hankel matrix. 
%
%

if nargin<5
    kcol= krow;
end

[n,s] = size(K);

H= zeros(krow*s,kcol*s);
Of = C;
Cp = K;

for j=1:krow 
    Of = [C;Of*A];
end;
for j=1:kcol
    Cp = [K,A*Cp];
end

H= Of*Cp;
