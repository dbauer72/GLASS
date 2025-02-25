s= 2;
n = 2;
Ts = [100,200,400,800,1600];
M = 1000; 

A = diag([0.7,.2]);
A0 = [A];
K0 = [eye(2)];
C0 = [A-eye(2)];

th0 = theta_urs();
th0.A = A0;
th0.K = K0;
th0.C = C0;
th0.D = zeros(s,0);
th0.B = zeros(n,0);
th0.Omega = eye(s);

im0 = impulse(th0);

thg = theta_urs();
thg.A= A;
thg.B = eye(s);
thg.C= A;
thg.D= eye(s);

h = waitbar(0,'Please wait'); 
parest_c = zeros(M,n^2+2*s*n,5);
parest_c2 = zeros(M,n^2+2*s*n,5);
parest_i = zeros(M,n^2+2*s*n,5);
im_norm = zeros(M,3,5);

im_c = zeros(M,12,5);
im_c2 = zeros(M,12,5);
im_i = zeros(M,12,5);

Abar_c = zeros(M,n^2,5);
Abar_c2 = zeros(M,n^2,5);
Abar_i = zeros(M,n^2,5);

for t=1:5
    T = Ts(t);
    for m=1:M
        waitbar(m/M);
        yg = idsim(T+1,thg);
        y = diff(yg);
        [result,thc,Ac,Kc,Cc,Omegac,thi,lle] = StSp_I0(y,s,n,10,.1);
        [result,thc2,Ac2,Kc2,Cc2,Omegac2,thi,lle] = StSp_I0(y,s,n,10,-1);
    
        % renormalize initial est
        [thi,Ai,BKi,Ci] = ss2ech_n(thi.A',thi.C',thi.K');
        thi.A= Ai';thi.K=Ci';thi.C = BKi';
        % parameter estimates 
        parest_c(m,:,t)= [Ac(:)',Kc(:)',Cc(:)'];
        parest_c2(m,:,t)= [Ac2(:)',Kc2(:)',Cc2(:)'];
        parest_i(m,:,t)= [thi.A(:)',thi.K(:)',thi.C(:)'];
    
        % matrices 
        Abarc = Ac-Kc*Cc;
        Abarc2 = Ac2-Kc2*Cc2;
        Abari = thi.A-thi.K*thi.C;
        
        Abar_c(m,:,t)= Abarc(:)';
        Abar_c2(m,:,t)= Abarc2(:)';
        Abar_i(m,:,t)= Abari(:)';

        % impulse responses 
        imc = impulse(thc);
        imc2 = impulse(thc2);
        imi = impulse(thi);

        for im=1:3
            im_c(m,(im-1)*4+[1:4],t)=reshape(squeeze(imc(:,:,im+1)),1,4);
            im_c2(m,(im-1)*4+[1:4],t)=reshape(squeeze(imc2(:,:,im+1)),1,4);
            im_i(m,(im-1)*4+[1:4],t)=reshape(squeeze(imi(:,:,im+1)),1,4);
        end
        im_norm(m,:,t)= [norm(imi-im0,'fro'),norm(imc-im0,'fro'),norm(imc2-im0,'fro')];
    end
end
close(h);


mi_no = squeeze(mean(im_norm,1));
save simu_over_diff;

% norm of impulse response error 
clf;
plot(Ts,(ones(3,1)*Ts).*(mi_no.^2))
set(gca,'fontsize',16)
xticks(Ts);
xtickangle(45)
xlabel('T')
ylabel('T*MNE')
legend({'CVA','qMLE','PEM'})

print -dpng MNE.png

% impulse response for T=200 
[xi,yi] = kdfft1(squeeze(im_i(:,1,2)),'knorm',1024,0.01);
[xc,yc] = kdfft1(squeeze(im_c(:,1,2)),'knorm',1024,0.01);
[xc2,yc2] = kdfft1(squeeze(im_c2(:,1,2)),'knorm',1024,0.01);

plot(xi,yi);
hold on;
plot(xc,yc);
plot(xc2,yc2);
set(gca,'fontsize',16)
xlabel('K_1(1,1)');
ylabel('kde')
xline(-.3,'Linewidth',2)
legend({'CVA','qMLE','PEM',''})

print -dpng kde_200.png

% impulse response for T=1600 
[xi,yi] = kdfft1(squeeze(im_i(:,1,5)),'knorm',1024,0.003);
[xc,yc] = kdfft1(squeeze(im_c(:,1,5)),'knorm',1024,0.003);
[xc2,yc2] = kdfft1(squeeze(im_c2(:,1,5)),'knorm',1024,0.003);

clf;
plot(xi,yi);
hold on;
plot(xc,yc);
plot(xc2,yc2);
set(gca,'fontsize',16)
xlabel('K_1(1,1)');
ylabel('kde')
xline(-.3,'Linewidth',2)
legend({'CVA','qMLE','PEM',''})

print -dpng kde_1600.png

% density for trace of Abar. 
% For T=200
j=2
[xi,yi] = kdfft1(squeeze(Abar_i(:,1,j))+squeeze(Abar_i(:,4,j)),'knorm',1024,0.01);
[xc,yc] = kdfft1(squeeze(Abar_c(:,1,j))+squeeze(Abar_c(:,4,j)),'knorm',1024,0.01);
[xc2,yc2] = kdfft1(squeeze(Abar_c2(:,1,j))+squeeze(Abar_c2(:,4,j)),'knorm',1024,0.01);

clf;
plot(xi,yi);
hold on;
plot(xc,yc);
plot(xc2,yc2);
set(gca,'fontsize',16)
xlabel('trace(Abar)');
ylabel('kde')
legend({'CVA','qMLE','PEM',''})

print -dpng abar_200.png

% for T=1600
trabar = squeeze(Abar_c2(:,1,j))+squeeze(Abar_c2(:,4,j));

j=5;
[xi,yi] = kdfft1(squeeze(Abar_i(:,1,j))+squeeze(Abar_i(:,4,j)),'knorm',1024,0.001);
[xc,yc] = kdfft1(squeeze(Abar_c(:,1,j))+squeeze(Abar_c(:,4,j)),'knorm',1024,0.001);
[xc2,yc2] = kdfft1(squeeze(Abar_c2(:,1,j))+squeeze(Abar_c2(:,4,j)),'knorm',1024,0.001);

clf;
plot(xi,yi);
hold on;
plot(xc,yc);
plot(xc2,yc2);
set(gca,'fontsize',16)
xlabel('trace(Abar)');
ylabel('kde')
legend({'CVA','qMLE','PEM',''})

print -dpng abar_1600.png


% calculate biased systems 


R = cal_cov(th0,500);

inc_p = zeros(100,12);
for m=1:100
    [the,Ap,Bp,Cp]=CVA_cov(R,n,2*m,2*m);
    inc_p(m,:)=[Ap(:)',Bp(:)',Cp(:)'];
end

% plot path for systems 
plot(2*[1:100],inc_p);
set(gca,'fontsize',16)
xlabel('f=p');

print -dpng syst_p.png

% plot deviations from limit time p/ p^2. 
inc_inf = [A0(:)',K0(:)',C0(:)'];
inc_p_dev = inc_p - ones(100,1)*inc_inf;
inc_p_dev(:,1:4)=inc_p_dev(:,1:4).*(2*([1:100]'*ones(1,4))).^2;
inc_p_dev(:,5:12)=inc_p_dev(:,5:12).*(2*([1:100]'*ones(1,8)));

plot(2*[1:100],inc_p_dev);
set(gca,'fontsize',16)
xlabel('f=p');

print -dpng d_syst_p.png



% sd_c = sqrt(var(parest_c));
% sd_c2 = sqrt(var(parest_c2));
% sd_i = sqrt(var(parest_i));
% 
% plot(sd_c,sd_i,'x')
% 
% % same for T=1000
% T=1000;
% M=100;
% h = waitbar(0,'Please wait'); 
% Mparest_c = zeros(M,n^2+2*s*n);
% Mparest_c2 = zeros(M,n^2+2*s*n);
% Mparest_i = zeros(M,n^2+2*s*n);
% Mim_norm = zeros(M,3);
% 
% for m=1:M
%     waitbar(m/M);
%     yg = idsim(T+1,thg);
%     y = diff(yg);
%     [result,thc,Ac,Kc,Cc,Omegac,thi,lle] = StSp_I0(y,s,n,10,1);
%     [result,thc2,Ac2,Kc2,Cc2,Omegac2,thi,lle] = StSp_I0(y,s,n,10,-1);
%     
%     % renormalize initial est
%     [thi,Ai,BKi,Ci] = ss2ech_n(thi.A',thi.C',thi.K');
%     thi.A= Ai';thi.K=Ci';thi.C = BKi';
%     % parameter estimates 
%     Mparest_c(m,:)= [Ac(:)',Kc(:)',Cc(:)'];
%     Mparest_c2(m,:)= [Ac2(:)',Kc2(:)',Cc2(:)'];
%     Mparest_i(m,:)= [thi.A(:)',thi.K(:)',thi.C(:)'];
%     
%     % impulse responses 
%     Mimc = impulse(thc);
%     Mimc2 = impulse(thc2);
%     Mimi = impulse(thi);
% 
%     Mim_norm(m,:)= [norm(Mimi-im0,'fro'),norm(Mimc-im0,'fro'),norm(Mimc2-im0,'fro')];
% 
% end
% close(h);
% 
% sd_c = sqrt(var(Mparest_c));
% sd_c2 = sqrt(var(Mparest_c2));
% sd_i = sqrt(var(Mparest_i));
% 
% plot(sd_c,sd_i,'x')
% 
% sum(sd_c)
% sum(sd_c2)
% sum(sd_i)
% 
