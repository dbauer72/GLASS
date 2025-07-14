function [th,par] = param2th_BLT(par,n,sf,index);
% converts parameters into th structures via system matrices for regional
% models.
% 
% SYNTAX: [th,par,dth] = param2th(par,n,sf,index); 
%
% INPUT: param ... dx1 parameter vector.
%        n     ... integer; state dimension.
%        sf    ... 2xinteger; dimension of regional vars + star vars
%        index  ... [ci,cist] common trends in syst and in star vars 
%
% OUTPUT: th ... theta structure corresponding to state space system.
%         par ... remaining parameters, if any.
%
% REMARK: Does not implement the derivative as of yet.  
%
% AUTHOR: dbauer, 2.8.2024.

si = sf(1);
sist = sf(2);
ci = index(1);
cist = index(2); 

% if nargout>2
%     np =length(par);
%     dA = zeros(n,n,np);
%     dK = zeros(n,s,np);
%     dC = zeros(s,n,np);
% end

% --- calc matrices for regional model ---
[thi,par,dth] = param2th(par,si,n,ci);

% --- D ---
parD = par(1:(si*sist));
par(1:(si*sist)) = [];

D = reshape(parD,si,sist);

% --- B ---
parB1 = par(1:(ci*(sist-cist)));
par(1:(ci*(sist-cist))) = [];

B1 = reshape(parB1,ci,sist-cist);

parB2 = par(1:((n-ci)*sist));
par(1:((n-ci)*sist)) = [];

B2 = reshape(parB2,n-ci,sist);

B = zeros(n,sist);
B(1:ci,cist+1:end)=B1;
B(ci+1:end,:)=B2; 

% --- Omega ---
Omega = fill_lowtri(par,si);

% --- fill into theta structure ---
th = theta_urs;
th.which = 'SS';
th.A = thi.A;
th.K=thi.K;
th.C=thi.C;
th.D = D;
th.B = B;
th.Omega = Omega; 

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
