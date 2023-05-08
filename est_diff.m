% script for testing state space estimation of MA(1) with differencing. 


M= 5000;
s=2; 

betaf = [-0.8:(-0.04):-1];

for b=1:length(betaf)
    beta = betaf(b)
    Tf = [100,500,1000,2000,3000,4000];
    for t=1:length(Tf)
        T = Tf(t)
        for m=1:M
            u= randn(T+1,s);
            y = u(2:end,:)+beta*u(1:end-1,:);
    
            nmax = floor(T^0.25)*4;
    
            k = aicest(y,s,nmax);
            k_all(b,t,m) = k;
            [the,ae,ke,ce,Omegae] = CCA(y,2,k,k,0);
    
            IMe = impulse(the,10);
            vIme(m,:) = reshape(IMe,1,10*s^2);
        end
        me(t,:) = mean(vIme);
        ve(t,:) = var(vIme);
    end
    
    true = [1,0,0,1,beta,0,0,beta,zeros(1,32)];;
    dev = me - ones(length(Tf),1)*true;
    DT = diag(sqrt(Tf));
    plot(dev'*DT)
    
    set(gca,'fontsize',16)
    legend(num2str(Tf'));
    
    % same for MSE
    rmse = sqrt(dev.^2 + ve.^2);
    figure
    plot(rmse'*DT)
    
    set(gca,'fontsize',16)
    legend(num2str(Tf'));
    
    scaled_rmse(:,b) = DT * sum(rmse,2);
end;

save diff_syst

