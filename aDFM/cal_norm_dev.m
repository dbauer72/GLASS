function norms= cal_norm_dev(th,th_chi,HT,lag);
% calculates the norm of the deviation between estimates in th and true
% system in th_chi.

norms = zeros(1,lag+1);

LD = th_chi.D*th_chi.D';
LDe = HT*th.D*th.Omega*th.D'*HT';

norms(1) = norm(LD-LDe);
for j=1:lag
    LD = th_chi.C*th_chi.A^(j-1)*th_chi.B*th_chi.B'*(th_chi.A^(j-1))'*th_chi.C';
    LDe = HT*th.C*th.A^(j-1)*th.B*th.Omega*th.B'*(th.A^(j-1))'*th.C'*HT';
    norms(1+j) = norm(LD-LDe);    
end
