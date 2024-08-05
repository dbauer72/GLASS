function [th,par] = param2th_RM(par,n,sf,index);
% converts parameters into th structures via system matrices for regional
% models.
% 
% SYNTAX: [th,par,dth] = param2th(par,s,n,c,restrict); 
%
% INPUT: param ... dx1 parameter vector.
%        n     ... integer; state dimension.
%        sf    ... integer; dimension of full system
%        index  ... vector of structural integers; see th2param_th for
%                   explanation.
%        restrict ... structure describing the restrictions. 
%
% OUTPUT: th ... theta structure corresponding to state space system.
%         par ... remaining parameters, if any.
%
% REMARK: Does not implement the derivative as of yet.  
%
% AUTHOR: dbauer, 2.8.2024.

si = index(1);
ci = index(2);
cist = index(3); 

% if nargout>2
%     np =length(par);
%     dA = zeros(n,n,np);
%     dK = zeros(n,s,np);
%     dC = zeros(s,n,np);
% end

% --- find number of params in C1 and K1 ---
nparc = (sf-ci)*ci; % p.l.t. form and orthonormality. 
parc = par(1:nparc);
par(1:nparc)=[];
C1 = par2ortho_plt(parc,sf,ci);
% if nargout>2 
%     for j=1:length(parc)
%        dC(:,1:c,j) = dpar2ortho_LR(parc,s,c,j);
%     end
% end

npark = si*ci;
% cur =nparc; 

% --- fill in parameters in K1 ---
park = par(1:npark);
par(1:npark)=[];

% --- first K_{i,c} ---
Kic = reshape(park(1:ci*si),ci,si);

% --- second K_{i,c}^* = Q_{i,c}^*R_{i,c}^* 
sist = sf-si;
if (ci>0)
    if (ci>cist)
        nparqst = ci*cist - cist*(cist+1)/2;
        parq = par(1:nparqst);
        par(1:nparqst)=[];

        Qist = par2ortho(parq,ci,cist);

        nparst = cist*sist - cist*(cist-1)/2;
        parpst = par(1:nparst);
        par(1:nparst)=[];
        Rist = zeros(cist,sf-si);

        for j=1:cist
            Rist(j,j:end)=parpst(1:(sist-j+1));
            park(1:(sist-j+1))=[];
        %     if nargout>2
        %         for jj=1:(s-j+1)
        %             dK(j,j+jj-1,cur+jj)= 1;
        %         end
        %         cur = cur+s-j+1;
        %     end
        end;
    else % then both must be equal
        Qist = eye(ci);
        nparst = ci*sist;
        parpst = par(1:nparst);
        par(1:nparst)=[];
        Rist = reshape(parpst,ci,sist);
    end
else
    Qist = eye(ci);
    Rist = zeros(ci,sist);
end




% ------ put parts together -----
K1 = [Kic,Qist*Rist]; 

cur = nparc+npark;

% --- remaining parameters for stationary part ---
parbull = par(1:(n-ci)*sf*2);
if nargout>2
    [Abull,Kbull,Cbull,dAKC]= param2mat_bull(parbull,sf,n-ci);
else
    [Abull,Kbull,Cbull]= param2mat_bull(parbull,sf,n-ci);
end
par(1:(n-ci)*sf*2)=[];

% --- fill in the submatrices ---
A = [eye(ci),zeros(ci,n-ci);zeros(n-ci,ci),Abull];
K = [K1;Kbull];
C = [C1,Cbull];

% --- fill into theta structure ---
th = theta_urs;
th.which = 'SS';
th.A = A;
th.K=K;
th.C=C;

% if nargout>2
%     % fill in 
%     for j=(cur+1):(np-length(par))
%         dA(:,:,j) = [zeros(c,n);zeros(n-c,c),squeeze(dAKC.A(:,:,j-cur))];
%         dC(:,:,j) = [zeros(s,c),squeeze(dAKC.C(:,:,j-cur))];
%         dK(:,:,j)= [zeros(c,s);squeeze(dAKC.K(:,:,j-cur))];
%     end  
%     
%     % write into dth
%     dth.A = dA;
%     dth.K = dK;
%     dth.C = dC; 
% end
