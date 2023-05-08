function [A,B,C,D,K,Omega,th,gr_syst] = param2syst_X(param,s,n,m,c,D0,restrict,Omega); 
% takes the parameter vector, the integers and the restrictions and
% converts it to system plus innovation variance.
%
% SYNTAX: [A,K,C,D,Omega,th,gr_syst] = param2syst(param,s,n,m,c,restrict); 
%
% INPUT: param ... dx1 parameter vector.
%        s     ... integer; dimension of endogenous vars.
%        n     ... integer; state dimension.
%        m     ... integer; dimension of exogenous vars.
%        c     ... integer; number of common trends.
%        restrict ... structure describing the restrictions. 
%
% OUTPUT: A,B,C,D,K ... state space system.
%         Omega   ... sxs innovation variance.
%         th      ... thetat structure corresponding to A,B,C,D,K.
%         gr_syst ... if wanted in the output contains the derivative of
%                       the system matrices with respect to the parameters. 
%                      Array: gr_syst(1).A, .K, .C, .D, .Omega contains the
%                      derivative of A,K,C,D, Omega respectively w.r.t. the
%                      first parameter, gr_syst(2).A, ... w.r.t. to the
%                      second param.
%
% REMARK: uses the parameterization corr. to the param paper. 
%         + only stationary case for now. 
%
% AUTHOR: dbauer, 27.4.2023

det_res = restrict.det_res;
if isfield(restrict,'scale')
    scale = diag(restrict.scale); 
else
    scale = eye(length(param));
end;

param = scale*param;


% extract parameters for Omega
if nargin<8 % Omega not contained
    parom = param(1:(s*(s+1)/2));
    nom = length(parom);
    param(1:(s*(s+1)/2))=[];    
else
    parom = extr_lowtri(Omega);
    nom = length(parom);
end

if nargout>6
    np = length(param)+nom;
    [Omega,dOmega] = fill_lowtri(parom,s);
else
    np = length(param);
    Omega = fill_lowtri(parom,s);
end

% convert params to matrices
nd = s*m;
if det_res 
    nd = nd- (s-c);
end;

if nargout>7 % not implemented yet!
    nth = length(param);
    [th,param,dth] = param2th_X(param,s,m,n,c,D0,restrict);
%    nth = nth-length(param);    
else
    [th,param] = param2th_X(param,s,m,n,c,D0,restrict);
end

A = th.A;
K = th.K;
C = th.C;
B = th.B;
D = th.D; 

% % parameters for deterministics. 
% if c>0
%     [Q,R]= qr(C(:,1:c));
% else
%     Q = eye(s);
% end;
% 
% if (det_res>0)&(c>0)
%     D = zeros(s,m);
%     pc = param(1:c);
%     D(:,1)=C(:,1:c)*pc(:);
% 
%     param(1:c)=[];
%     if m>1
%         D(:,2:end)= reshape(param,s,m-1);
%     end
%     if nargout>6
%         for j=1:c
%             dD(:,1,j)=C(:,j);
%         end
%         if m>1 
%             dp = eye(s*(m-1));
%             for j=1:(s*(m-1))
%                 dD(:,:,j) = reshape(dp(:,j),s,m-1);
%             end
%         end
%     end
% else
%     D = reshape(param,s,m);
%     if nargout>6
%         dp = eye(s*m);
%         for j=1:s*m,
%             dD(:,:,j) = reshape(dp(:,j),s,m);
%         end
%     end;
% end

th.Omega =Omega;
th.ur = 'I(0)';

if nargout>7 % gradient also wanted: not implemented yet.   
    for j=1:np
        da = zeros(n,n);
        dk = zeros(n,s);
        db = zeros(n,m);
        dc = zeros(s,n);
        dd = zeros(s,m);
        dom = zeros(s,s);
        
        if j<= nom
            dom = squeeze(dOmega(:,:,j));
        end
        if (j> nom)&(j<=nom+nth) % parameter for (A,B,C,D,K)
            da = squeeze(dth.A(:,:,j-nom));
            db = squeeze(dth.B(:,:,j-nom));
            dc = squeeze(dth.C(:,:,j-nom));
            dd = squeeze(dth.D(:,:,j-nom));
            dk = squeeze(dth.K(:,:,j-nom));
        end
        
        gr_syst(j).A= da;    
        gr_syst(j).B= db;
        gr_syst(j).C= dc;    
        gr_syst(j).K= dk;    
        gr_syst(j).D= dd;   
        gr_syst(j).Omega= dom;    
    end
end


