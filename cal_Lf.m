function [Lf,S]= cal_Lf(th,f);
% cal_Lf calculates the matrix Lf of Anderson and Deistler with f block
% rows.
%
% SYNTAX:  [Lf]= cal_Lf(th,f);
% 
% INPUT:  th ... theta structure
%         f  ... integer; number of blocks. 
% OUTPUT: Lf ... matrix.
%
% AUTHOR: dbauer, 14.6.2024. 

[r,q]= size(th.D); 
n = size(th.A,1);

Lf = zeros(f*r,f*q+n);

A = th.A;B= th.B;C=th.C;D=th.D; 

for j= 1:f
    Lf((j-1)*r+[1:r],1:n) = C*A^(j-1); 
    Lf((j-1)*r+[1:r],(j-1)*q+[1:q]+n)= D;
    if (j>1)
        for p = 1:(j-1)
            Lf((j-1)*r+[1:r],(p-1)*q+[1:q]+n) = C*A^(j-p-1)*B;
        end
    end
end

[~,s,~] = svd(Lf);
S = diag(s);

