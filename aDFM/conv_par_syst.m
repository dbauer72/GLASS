function [Ai,Bi,Ci,Di,Psi]= conv_par_syst(theta,s,n);
%
%

Ai= reshape(theta(1:n^2),n,n);
theta(1:n^2)=[];
Bi = reshape(theta(1:(n*s)),n,s);
theta(1:(n*s))=[];
Ci = reshape(theta(1:(n*s)),s,n);
theta(1:(n*s))=[];

Di= zeros(s,s);
for j=1:s
    Di(j,1:j)=theta(1:j);
    theta(1:j)=[];
end;

Psi = fill_lowtri(theta,s);
