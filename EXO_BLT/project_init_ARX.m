function th = project_init_ARX(thi,si,sist,ni,nist,ci,cist);
% project_init_ARX projects an estimated ARX model onto 
% the GLASS regional model.
% 
% SYNTAX: th = project_init_ARX(thi,si,sist,ci,cist);
% 
% INPUTS:  thi ... initial estimate of ARX system
%          si  ... integer; dimension of y_{i,t}.
%          sist ... integer; dimension of star var y_{i,t}^*
%          ci  ... integer; number of common trends in y_{i,t}
%          cist ... integer; number of common trends in y_{i,t}^*.
%
%  
% OUTPUT: th ... theta structure 
%
% AUTHOR: dbauer, 14.7.2025.

% extract p and q.
a = thi.a; 

phi = a(:,(si+1):end);
p = floor(size(phi,2)/si+0.1);

psif = thi.d; 
D= psif(:,1:sist);
psi = psif(:,(sist+1):end);

q = floor(size(psi,2)/sist+0.1);
pq = q; 

if p<q
    pq = q;
    phi(:,p*si+[1:(q-p)*si])=0;
end
if q<p 
    pq = q; 
    psi(:,q*si+ [1:(p-q)*sist])=0;
end

% calculate companion matrix for a(z)
Ky = size(phi,2); 
Kx = size(psi,2);
K = Ky+Kx; 

Pi = zeros(K,K);
Pi(1:si,:)=[-phi,psi];
Pi((si+1):Ky,1:(Ky-si))=eye(Ky-si);
Pi((Ky+sist+1):end,(Ky+1):(end-sist))=eye(Kx-sist);

Bi = zeros(K,sist);
Bi(Ky+[1:sist],:)=eye(sist);
Ki = zeros(K,si);
Ki(1:si,1:si)=eye(si);

Ci = Pi(1:si,:);

% calculate Hankel matrix
H = my_hank(Pi,[Bi,Ki],Ci,5*ni,5*ni);

% different rewriting 
% Ai2 = zeros(Ky,Ky);
% Bi2 = zeros(Ky,sist);
% Ki2 = zeros(Ky,si);
% Ci2 = zeros(si,Ky);
% 
% Ci2(:,1:si)=eye(si);
% Ai2(1:(Ky-si),(si+1):end)=eye(Ky-si);
% for j=1:pq
%     Ai2((j-1)*si+[1:si],1:si) = -phi(:,(j-1)*si+[1:si]);
%     Bi2((j-1)*si+[1:si],1:sist) = psi(:,(j-1)*sist+[1:sist])-phi(:,(j-1)*si+[1:si])*D;
%     Ki2((j-1)*si+[1:si],1:si) = -phi(:,(j-1)*si+[1:si]);    
% end
% 
% 
% IR = impulse(thi,10*ni);



% approximate Hankel matrix
[u,s,v]= svd(H);
plot(diag(s),'x');

Of = u(:,1:ni);
Cp = Of'*H;

% extract system matrices and normalize system
[Jit,Bit,Cit,Kit] = calculate_normed_system_from_OfCp(Of,Cp,D,si,sist,ni,nist,ci,cist);

% store results 
th = theta_urs();
th.which ='SS';
th.A = Jit;
th.B = Bit;
th.C = Cit;
th.D = D;
th.K = Kit;
th.Omega = thi.Omega;

