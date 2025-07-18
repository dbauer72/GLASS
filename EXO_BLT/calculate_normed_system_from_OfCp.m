function [Jit,Bit,Cit,Kit] = calculate_normed_system_from_OfCp(Of,Cp,D,si,sist,ni,nist,ci,cist);
%  calculate_normed_system_from_OfCp.
%
% SYNTAX: [A,B,C,K] = calculate_normed_system_from_OfCp(Of,Cp,D,si,sist,ni,nist,ci,cist);
%
% INPUTS: Of ... observability matrix
%         Cp ... controllability matrix.
%         D  ... matrix.
%         si, sist ... integers, dimensions of y and y*
%         ni,nist ... integers, dimensions of states for y and y*
%         ci, cist... integers, number of common trends for y and y*
%
%  OUTPUT: (A,B,C,K) ... system in Bauer Wagner canonical form.
%
% AUTHOR: dbauer, 14.7.2025. 


Ci = Of(1:si,1:ni);
Ki = Cp(:,(sist+1):(sist+si)); 
Bi = Cp(:,1:sist) + Ki*D;

Ai = Of(1:(end-si),:)\Of((si+1):end,:);

% project to get number of common trends right. 
[v,d]= eig(Ai);
ds = abs(diag(d));
[dsort,I] = sort(ds);
[dsort,I]=sort(ds,'descend');

% calculate transformation and normalisation 
Trafo = v(:,I);
Jit = inv(Trafo)*Ai*Trafo;

% unit eigenvalue part
Jit(1:ci,1:ci)=eye(ci);
Jit(1:ci,(ci+1):end) = 0;
Jit((ci+1):end,1:ci) = 0;
Cit = Ci*Trafo;
Bit = inv(Trafo)*Bi;
Kit = inv(Trafo)*Ki;

Cit(:,1:ci)=real(Cit(:,1:ci));
Bit(1:ci,:)=real(Bit(1:ci,:));
Kit(1:ci,:)=real(Kit(1:ci,:));

[Q,R]= qr(Cit(:,1:ci));
Cit(:,1:ci)=Q(:,1:ci);
Bit(1:ci,:)=R(1:ci,1:ci)*Bit(1:ci,:);
Kit(1:ci,:)=R(1:ci,1:ci)*Kit(1:ci,:);

% stable part
ind = (ci+1):ni;

% extract matrices
Abull = Jit(ind,ind);
Cbull = Cit(:,ind);
Kbull = [Bit(ind,:),Kit(ind,:)];

% calculate impulse response 
IR =zeros(si,sist+si,2*ni);
CA = Cbull;
for j=1:(2*ni)
    IR(:,:,j)=real(CA*Kbull);
    CA = CA*Abull;
end

% calculate Hankel matrix
H = zeros(si*ni,(ni+1)*(si+sist));
for a=1:ni
    for b=1:(ni+1)
        H((a-1)*si+[1:si],(b-1)*(sist+si)+[1:(sist+si)])= squeeze(IR(:,:,a+b-1));
    end
end

[u,s,v]= svd(H);
nbull= ni- ci;
Of = u(:,1:nbull);
Cp = Of'*H;
Cbull = Of(1:si,:);
Abull = Of(1:(end-si),:)\Of((si+1):end,:);
Bbull = Cp(:,1:sist);
Kbull = Cp(:,sist+[1:si]);

% fill in the matrices 
Jit(ind,ind)= Abull;
Cit(:,ind) = Cbull;
Bit(ind,:)=Bbull;
Kit(ind,:) = Kbull;

% restriction to I(1)
Bit(1:ci,1:cist)=0;

%Bit = Bit + Kit*D;