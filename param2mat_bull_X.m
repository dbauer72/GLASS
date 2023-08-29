function [Abull,Bbull,Cbull,Kbull,dABCK]= param2mat_bull_X(parbull,s,m,n,D0);
% writes the parameters into the matrices for system in echelon form,
% including exogenous input matrix B. 
%
% SYNTAX: [Abull,Bbull,Cbull,Kbull,dABCK]= param2mat_bull_X(parbull,s,m,n,D0);
%
% INPUT: parbull ... dx1 parameter vector.
%        s       ... integer; output dimension.
%        m       ... integer; dimension exogenous vars.
%        n       ... integer; state dimension.
%        D0  ... indicator; If D0=1 (TRUE), then D=0 is used. 
%
% OUTPUT: Abull,Bbull,Cbull, Kbull ... system matrices.
%         dABCK   ... if nargout>3, then also derivatives of the system
%         matrices with respect to the parameters are calculated. 
%                     dABCK.A: array of dimension nxnxd. 
%                     dABCK.K: array of dimension nxsxd. 
%                     dABCK.C: array of dimension sxnxd. 
%                     dABCK.B: array of dimension nxmxd.
%
% REMARK: + only stationary systems up to now. 
%
% AUTHOR: dbauer, 27.4.2023

if nargout>3
    np = length(parbull);
    dA = zeros(n,n,np);
    dK = zeros(n,s,np);
    dC = zeros(s,n,np);
    dB = zeros(n,m,np);    
end

% if s>=n -> parameters for all matrices contained. 
if s>=n
    % parameters for C: (s-n)*n
    parc= parbull(1:(s-n)*n);
    Cbull = [eye(n);reshape(parc,s-n,n)];
    parbull(1:(s-n)*n)=[];
    % parameters for A:
    Abull = reshape(parbull(1:n^2),n,n);
    parbull(1:n^2)=[];
    
    % parameters for Kbull
    Bbull = reshape(parbull(1:n*m),n,m);
    parbull(1:n*m)=[];
    Kbull = reshape(parbull(1:n*s),n,s);

    if nargout>3
        dcp = eye((s-n)*n);
        for j=1:(s-n)*n
            dC(n+[1:(s-n)],:,j)=reshape(dcp(:,j),s-n,n);
        end;
        dap = eye(n^2);
        for j=1:n^2
            dA(:,:,j+(s-n)*n)=reshape(dap(:,j),n,n);
        end;
        dbp = eye(n*m);
        for j=1:n*m
            dB(:,:,s*n+j)=reshape(dbp(:,j),n,m);
        end                   
        dkp = eye(n*s);
        for j=1:n*s
            dK(:,:,s*n+n*m+j)=reshape(dkp(:,j),n,s);
        end                   
    end
else
    % if s<n -> no parameter is C. 
    Cbull = [eye(s),zeros(s,n-s)];
    
    % A contains s*n parameters. 
    Abull = zeros(n,n);
    Abull(1:(n-s),s+1:n)=eye(n-s);
    para = parbull(1:n*s);
    Abull(n-s+1:end,:)=reshape(para,s,n);
    
    parbull(1:n*s)=[];
    
    % B,K contains the remaining parameters. 
    % parameters for Kbull
    Bbull = reshape(parbull(1:n*m),n,m);
    parbull(1:n*m)=[];
    Kbull = reshape(parbull(1:n*s),n,s);

    if nargout>3
        dap = eye(n*s);
        for j=1:n*s
            dA(n-s+1:end,:,j)=reshape(dap(:,j),s,n);
        end;
        dbp = eye(n*m);
        for j=1:n*m
            dB(:,:,s*n+j)=reshape(dbp(:,j),n,m);
        end                   
        dkp = eye(n*s);
        for j=1:n*s
            dK(:,:,s*n+n*m+j)=reshape(dkp(:,j),n,s);
        end                 
    end
end;

% % switch roles
% Abull = Abull';
% Cb = Cbull;
% Cbull = Kbull';
% Kbull = Cb'; 
% 
if nargout>3
    dABCK.A = dA; %permute(dA,[2,1,3]);
    dABCK.B = dB; 
    dABCK.K = dK; %permute(dC,[2,1,3]);
    dABCK.C = dC; %permute(dK,[2,1,3]); 
end
