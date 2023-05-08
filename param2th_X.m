function [th,par,dth] = param2th_X(par,s,m,n,c,D0,restrict);
% converts parameters into th structures via system matrices.
% same as param2th but without exogenous inputs
% 
% SYNTAX: [th,par,dth] = param2th_X(par,s,m,n,c,D0,restrict);
%
% INPUT: param ... dx1 parameter vector.
%        s     ... integer; dimension of endogenous vars.
%        m     ... integer; dimension of exogenous vars.
%        n     ... integer; state dimension.
%        c     ... integer; number of common trends.
%        D0  ... indicator; If D0=1 (TRUE), then D=0 is used. 
%        restrict ... structure describing the restrictions. 
%
% OUTPUT: th ... theta structure corresponding to state space system.
%         par ... remaining parameters for D. 
%         dth      ... if wanted in the output contains the derivative of
%                       the system matrices with respect to the parameters. 
%                      Array: dth.A: array nxnxd of derivatives of A w.r.t.
%                      the parameters. 
%
% REMARK: uses the parameterization corr. to the param paper. 
%         + only stationary case for now. 
%
% AUTHOR: dbauer, 27.4.2023

if nargin<7
    restrict = [];
end

if nargin<6
    D0=0; % default: D0 is estimated. 
end


if nargout>2
    np =length(par);
    dA = zeros(n,n,np);
    dK = zeros(n,s,np);
    dC = zeros(s,n,np);
    dD = zeros(s,m,np);
    dB = zeros(n,m,np);
end

% I(1) not implemented yet. 
% % --- find number of params in C1 and K1 ---
% if isfield(restrict,'H')
%     nparc =restrict.nparc;
%     parc = par(1:nparc);
%     par(1:nparc)=[];
%     C1 = feval(restrict.par2ortho,parc,s,c,restrict.H);
% else
%     nparc = s*c-c*(c+1)/2;
%     parc = par(1:nparc);
%     par(1:nparc)=[];
%     C1 = par2ortho_LR(parc,s,c);
%     if nargout>2 
%         for j=1:length(parc)
%             dC(:,1:c,j) = dpar2ortho_LR(parc,s,c,j);
%         end
%     end
% end
% 
% npark = s*c - c*(c-1)/2;
% cur =nparc; 
% 
% % --- fill in parameters in K1 ---
% park = par(1:npark);
% K1 = zeros(c,s);
% for j=1:c
%     K1(j,j:end)=park(1:(s-j+1));
%     park(1:(s-j+1))=[];
%     if nargout>2
%         for jj=1:(s-j+1)
%             dK(j,j+jj-1,cur+jj)= 1;
%         end
%         cur = cur+s-j+1;
%     end
% end;
% 
% par(1:npark)=[];
% cur = nparc+npark;

% --- remaining parameters for stationary part ---
np = 2*n*s + n*m; 
parbull = par(1:np);
if nargout>2
    [Abull,Bbull,Cbull,Kbull,dABCK]= param2mat_bull_X(parbull,s,m,n-c,D0);
else
    [Abull,Bbull,Cbull,Kbull] = param2mat_bull_X(parbull,s,m,n-c,D0);
end
par(1:np)=[];

% --- fill in the submatrices ---
% A = [eye(c),zeros(c,n-c);zeros(n-c,c),Abull];
% K = [K1;Kbull];
% C = [C1,Cbull];

A= Abull;
C= Cbull;
B = Bbull;
K = Kbull;

% --- handle D ---
D = zeros(s,m); 
if (D0 ~= 1)
    parD = par(1:s*m);
    D = reshape(parD,s,m);
end

% --- fill into theta structure ---
th = theta_urs;
th.which = 'SS';
th.A = A;
th.K=K;
th.C=C;
th.B=B;
th.D=D;


cur = 0; 
if nargout>2
    % fill in 
    for j=(cur+1):(np-length(par))
        dA(:,:,j) = [zeros(c,n);zeros(n-c,c),squeeze(dABCK.A(:,:,j-cur))];
        dC(:,:,j) = [zeros(s,c),squeeze(dABCK.C(:,:,j-cur))];
        dK(:,:,j)= [zeros(c,s);squeeze(dABCK.K(:,:,j-cur))];
        dB(:,:,j) = [squeeze(dABCK.B(:,:,j-cur))];
    end  
    cur = np-length(par);
    for j=1:length(par)
        ddp = eye(s*m);
        for j=1:s*m
            dD(:,:,cur+j)=reshape(ddp(:,j),s,m);
        end   
    end

    % write into dth
    dth.A = dA;
    dth.B = dB;
    dth.K = dK;
    dth.C = dC; 
    dth.D = dD; 
end
