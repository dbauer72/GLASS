% skript for estimating global State Space Models. 
% D. Bauer, 26.4.2023. 


% construct the system: 4 regions, one variable in each.  

% weights for star variables include an average of adkacent regions: 
Wi = zeros(4,4);
Wi(1,2)=1;
Wi(2,[1,3])=0.5;
Wi(3,[2,4])=0.5;
Wi(4,3)=1;

% common variables: average over all.
Wom = ones(1,4); 

% diagonal system: Ci=Ki=1, bAi = 0.1+0.1*i; 
C= eye(4);
bA = diag([0.2,0.3,.4,.5]); 

K = (eye(4)+ diag(randn(4,1))*Wi + randn(4,1)*Wom*0)*0.25;

while (any(abs(eig(bA+K*C))==true)) 
    K = (eye(4)+ diag(randn(4,1))*Wi + randn(4,1)*Wom*0)*0.25;
end

A = bA+K*C;
% true values for B and D=0: 
for j=1:4
    Xi = zeros(3,4);
    Xi(1,j)=1;
    Xi(2,:)=Wi(j,:);
    Xi(3,:)=Wom;

    theta(j,:) = K(j,:)*Xi'*inv(Xi*Xi');
end

% D=0? -> D0=1 (TRUE). 
D0=1;


% simulate the process
T = 5000;
th = ss2ech_n(A,K,C);
u = randn(4,T);
x = zeros(4,1);
y = zeros(4,T);
for t=1:T
    y(:,t) = C*x + u(:,t);
    x = A*x + K*u(:,t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% estimate with true Wi   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% cycle over regions to estimate the model for true star and common
% variables
nj = 1;


for j=1:4
    ystar = Wi(j,:)*y;
    ycom = Wom*y;

    z = [y(j,:)',ystar'];
    [thi(j),Ai,Bi,Ci,Di,Ki,Omegai] = CCAX(z,1,nj,10,10,1,D0);
    [result(j),thc(j),Ac,Bc,Cc,Dc,Kc,Omegac,~,lle] = StSp_I0X(z,1,1,10,D0,-1);
end

res = u'*0;
for j=1:4
    res(:,j)= result(j).res;
end

% A0 
A0 = zeros(4,4);
n=0;
for j=1:4
    n = n+size(thc(j).A,1);
    D = thc(j).D;
    A0(j,:)= D(1)*Wi(j,:); %+D(2)*Wom;
end

Af = zeros(n,n);
Kf = zeros(n,4);
Cf = Kf';

n= 0;
for j=1:4
    nj = size(thc(j).A,1);
    % bar A: inverse dynamic matrix
    Af(n+[1:nj],n+[1:nj])=thc(j).A - thc(j).K*thc(j).C;
    % C: block diag
    Cf(j,n+[1:nj])=thc(j).C;
    % Kf: first part block diag. 
    Ki = thc(j).K;
    Kf(n+[1:nj],j)=Ki;
    % second part due to star vars 
    Kf(n+[1:nj],:)= Kf(n+[1:nj],:) + (thc(j).B(:,1)-Ki*thc(j).D(:,1))*Wi(j,:);
    % third part due to common vars. 
    %Kf(n+[1:nj],:)= Kf(n+[1:nj],:) + (thc(j).B(:,2)-Ki*thc(j).D(:,2))*Wom;
    
    % increment block number. 
    n = n+nj; 
end

Af = Af+Kf*Cf;
Cf = inv(eye(4)-A0)*Cf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  new trial with different Wi  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nj = 1;
% D=0? -> D0=1 (TRUE). 


Wi = ones(4,4)-eye(4);
for j=1:4
    ystar = Wi(j,:)*y;
    ycom = Wom*y;

    z = [y(j,:)',ystar'];
    [thi(j),Ai,Bi,Ci,Di,Ki,Omegai] = CCAX(z,1,nj,10,10,1,D0);
    [result1(j),thc1(j),Ac,Bc,Cc,Dc,Kc,Omegac,~,lle] = StSp_I0X(z,1,1,10,D0,-1);
end

res1 = u'*0;
for j=1:4
    res1(:,j)= result1(j).res;
end

% A0 
A01 = zeros(4,4);
n=0;
for j=1:4
    n = n+size(thc1(j).A,1);
    D = thc1(j).D;
    A01(j,:)= D(1)*Wi(j,:); %+D(2)*Wom;
end

Af1 = zeros(n,n);
Kf1 = zeros(n,4);
Cf1 = Kf1';

n= 0;
for j=1:4
    nj = size(thc1(j).A,1);
    % bar A: inverse dynamic matrix
    Af1(n+[1:nj],n+[1:nj])=thc1(j).A - thc1(j).K*thc1(j).C;
    % C: block diag
    Cf1(j,n+[1:nj])=thc(j).C;
    % Kf: first part block diag. 
    Ki = thc1(j).K;
    Kf1(n+[1:nj],j)=Ki;
    % second part due to star vars 
    Kf1(n+[1:nj],:)= Kf1(n+[1:nj],:) + (thc1(j).B(:,1)-Ki*thc1(j).D(:,1))*Wi(j,:);
    % third part due to common vars. 
    %Kf(n+[1:nj],:)= Kf(n+[1:nj],:) + (thc(j).B(:,2)-Ki*thc(j).D(:,2))*Wom;
    
    % increment block number. 
    n = n+nj; 
end

Af1 = Af1+Kf1*Cf1;
Cf1 = inv(eye(4)-A01)*Cf1;


% compare the outcome.
nd = zeros(2,20);
for j=1:20
    nd(1,j) = norm(Cf*Af^(j-1)*Kf-C*A^(j-1)*K);
    nd(2,j) = norm(Cf1*Af1^(j-1)*Kf1-C*A^(j-1)*K);
end

figure
plot(nd');

