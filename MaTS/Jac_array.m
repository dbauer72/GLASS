function [vPsi,A,B,dA,dB,dvPsi,JA,JB,JvPsi] = Jac_array(param,M,N,L,p);
%  calculates the Jacobian of the term evaluation 

[A,B,dA,dB] = param2MaTS(param,M,N,L,p);

vPsi = zeros(M*N,M*N,L); 
for l=1:L
    for jr=1:p
        vPsi(:,:,l) = vPsi(:,:,l) + kron(squeeze(B(:,:,l,jr)),squeeze(A(:,:,l,jr)));
    end
end

dvPsi = zeros(M*N,M*N,L,length(param));
for d=1:length(param)
    for l=1:L
        for jr=1:p
            dvPsi(:,:,l,d) = dvPsi(:,:,l,d) + kron(squeeze(B(:,:,l,jr)),squeeze(dA(:,:,l,jr,d))) + kron(squeeze(dB(:,:,l,jr,d)),squeeze(A(:,:,l,jr)));
        end
    end
end


JA = dA*0;
JB = dB*0;
JvPSi = dvPsi*0;

% now do the pertubation in each entry 
eps = 0.00001;
for j=1:length(param)
    param_n = param;
    param_n(j) = param(j) + eps;

    [An,Bn] = param2MaTS(param_n,M,N,L,p);
    JA(:,:,:,:,j) = (An-A)/eps;
    JB(:,:,:,:,j) = (Bn-B)/eps;

    vPsin = zeros(M*N,M*N,L);
    for l=1:L
        for jr=1:p
            vPsin(:,:,l) = vPsin(:,:,l) + kron(squeeze(Bn(:,:,l,jr)),squeeze(An(:,:,l,jr)));
        end
    end
    JvPsi(:,:,:,j) = (vPsin-vPsi)/eps;
end

