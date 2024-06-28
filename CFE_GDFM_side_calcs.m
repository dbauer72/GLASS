M= 100;
for m=1:M
    ths(1) = theta_urs();
    ths(1).A = 1.4*(rand(1)-.5);
    ths(1).B = 1;
    ths(1).C = ths(1).A;
    ths(1).D = 1;
    ths(1).Omega = .1;
    for j=2:N
        ths(j) = theta_urs();
        ths(j).A = 1.4*(rand(1)-.5);
        ths(j).B = 1;
        ths(j).C = ths(j).A;
        ths(j).D=1;
        ths(j).Omega = .5;
    end

    % common factor part.
    n = 3;
    r = 4;
    q = 2;
    p=2;
    th_chi = theta_urs();
    th_chi.C = randn(r,n);
    th_chi.A = diag(0.8*(rand(n,1)-.5));
    th_chi.A = diag([-.8,.4,.8]);
    %dA = diag(th_chi.A);
    %dA(abs(dA)<0.5)=0.5;
    %th_chi.A = diag(dA);
    th_chi.B = randn(n,q);
    th_chi.D = randn(r,q);
    th_chi.Omega = eye(q);

    [y,chi,u,x,e,F]= simu_GDFM(T,ths,th_chi,Lambda);
    %[Ft,rest,Lambdahat] = PCA_select_r(y,r);
    [th,A,K,C,D,Omega,nhat,se,x] = CCA_sing(F,3,2,2,2,0);

    svs(m,:)=se(:)';
end
