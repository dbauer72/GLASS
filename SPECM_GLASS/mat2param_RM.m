function param1 = mat2param_RM(C1,K1,index);
% converts C1 and K1 matrices into parameters assuming in canonical form
% specified by the index = [s(i),c(i),c(i)^*]:
%   s(i): dimension of y_{i,t}
%   c(i): number of common trends in i-th regional model.
%   c(i)^*: rank of K_{i,c}^* indicating the number of common trends spilling over.
% 
%
% SYNTAX: param1 = mat2param(C1,K1);
%
% INPUT: C1 ... sxc matrix;
%        K1 ... cxs matrix;
%
% OUTPUT: param ... dx1 parameter vector.
%
% REMARK: assumes that C1 and K1 correspond to the generic neighborhood of
% the canonical form such that C1*C1 = I_c and p.l.t. and K1
% =[K_{i,c},K_{i,c}^*] where K_{i,c}* = Q_{i,c}^* R_{i,c}^*.
%
% AUTHOR: dbauer, 2.8.2024.

si = index(1);
ci = index(2);
cist = index(3); 

% collect parameters in C_{i,c}
pc = ortho2par_plt(C1);

% collect parameters in K_{i,c} 
Ki = K1(:,1:si);
pk = Ki(:); 

% collect parameters in K_{i,c}^* 
pkst = [];
if (cist>0) 
    if (ci>cist) % K_{i,st} is non-zero and restricted. 
        [Qi,Ri]=qr(K1(:,si+1:end));

        Qist = Qi(:,1:cist);
        Rist = Ri(1:cist,:);

        pq = ortho2par(Qist);
        pr = [];
        for j=1:cist
            pr = [pr,Rist(j,j:end)];
        end
        pkst = [pq(:);pr(:)];
    else % not rank restricted 
        Qist = eye(ci);
        Rist = K1(:,si+1:end);
        pkst = Rist(:);
    end
end

param1 = [pc(:);pk(:);pkst];
