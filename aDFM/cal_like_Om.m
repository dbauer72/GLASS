function crit = cal_like_Om(Omegat); 
%
%

T = size(Omegat,3);
r = size(Omegat,1);

crit = 0;
for t=1:T
    crit = crit + log(det(squeeze(Omegat(:,:,t))));
end