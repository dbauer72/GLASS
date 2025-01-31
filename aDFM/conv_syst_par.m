function theta = conv_syst_par(Ai,Bi,Ci,Di,Psi);
%
%

[n,s]=size(Bi);

vDi=[];
for j=1:s
    vDi = [vDi;Di(j,1:j)'];
end

thPsi = extr_lowtri(Psi);
theta= [Ai(:);Bi(:);Ci(:);vDi(:);thPsi(:)];
