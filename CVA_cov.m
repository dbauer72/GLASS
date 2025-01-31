%function [th,pos]=racca(y,n,k)
%
% realization-based algorithm  
% cooordinate system of CCA
% pos = 1, if procedure failed; else pos=0

function [th,Ap,Bp,Cp]=CVA_cov(R,n,f,p);

s = size(R,1);

Gp=zeros(s*p,s*p);
Gf=zeros(s*f,s*f);
Hfp=zeros(s*f,s*p);

for i=1:f
    for j=1:p
        Hfp((i-1)*s+[1:s],(j-1)*s+[1:s])=squeeze(R(:,:,i+j));
    end
end

for i=1:f
    for j=1:i
        Gf((i-1)*s+[1:s],(j-1)*s+[1:s])=squeeze(R(:,:,i-j+1));
        Gf((j-1)*s+[1:s],(i-1)*s+[1:s])=squeeze(R(:,:,i-j+1))';
    end    
end

for i=1:p
    for j=1:i
        Gp((i-1)*s+[1:s],(j-1)*s+[1:s])=squeeze(R(:,:,i-j+1))';
        Gp((j-1)*s+[1:s],(i-1)*s+[1:s])=squeeze(R(:,:,i-j+1));
    end    
end

%  ---------------  coordinate system of CCA   ------------------
W1=chol(Gf);W2=chol(Gp);
iW1=inv(W1);iW2=inv(W2);
[U,S,V]=svd(iW1'*Hfp*iW2);
%s1=diag(S);s1=s1(1:n);sq=diag(1 ./sqrt(s1));

Kp = V(:,1:n)'*inv(W2)';
Pp = Kp*Gp*Kp';
Ap=Kp*[Hfp(1:s,:);Gp(1:(end-s),:)]*Kp'*inv(Pp);
Cp = Hfp(1:s,:)*Kp'*inv(Pp); 
Omegap= squeeze(R(:,:,1))- Cp*Pp*Cp'; 
M = Kp*[squeeze(R(:,:,1)),Hfp(1:s,1:(end-s))]';
Bp=(M - Ap*Pp*Cp')*inv(Omegap); 


% ----------   transformation to Echelon form
% generic neighbourhood
th = theta_urs();
th.A=  Ap;
th.K = Bp;
th.C = Cp;
th.Omega = Omegap; 

[th,Ap,Bp,Cp] = ss2ech_n(th.A',th.C',th.K');
th.A= Ap';th.K=Cp';th.C = Bp';
Ap = th.A;
Bp = th.K;
Cp = th.C; 
