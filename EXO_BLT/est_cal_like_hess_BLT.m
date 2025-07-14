function [pare,llc2,resc2] = est_cal_like_hess_BLT(z,n,sf,index,thi,restrict);
% est_cal_like_hess_BLT performs minimization of the quasi likelihood taking
% illconditioning of the Hessian into account. 
% Starting values are contained in thi.
%
% SYNTAX: pare = est_cal_like_hess_BLT(z,n,sf,index,thi,restrict);
%
% INPUTS:  z ... Tx(s+m); observations
%        n    ... integer; state dimension
%        sf   ... pair of integers; dimension vars region i; dim ystar. 
%        index    ... vector of integers; see 'cal_quasi_like_RM' for
%                   detailed explanation
%          thi ... theta structure containing the initial estimates.
%          restrict ... restrict structure used for calculating the system
%                  for given parameters. 
% 
% OUTPUTS:  pare ... dx1 optimal parameter value
%
% REMARKS: + optimization is performed repeatedly (5 iterations hard
% coded!) where in each iteration the parameters are reweighted by inverses
% of the square roots of the diagonal values of the estimated Hessian. This
% makes the problem more well conditioned. 
%         + values of rescaling are passed on to the parameterization in
%         restrict.scale. 
% 
% AUTHOR: dbauer, 14.1.2020.
parami = th2param_BLT(thi,index,1);

% improve estimate
%options = optimoptions('fminunc','display','iter');
options = optimoptions('fminunc','display','iter','MaxIterations',1000);
options.MaxFunctionEvaluations = 20000;
%restrict.scale = ones(length(parami),1);
%restrict.Omega = thi.Omega; 

Pbull = -1; 

[pare,fval,exitflag] = fminunc(@(x) cal_quasi_like_BLT(x,z,n,sf,index,Pbull,restrict),parami,options);
[llc2,resc2] = cal_quasi_like_BLT(pare,z,n,sf,index,Pbull,restrict);

for j=1:5
    
    %Omegai = resc2'*resc2/size(resc2,1);
    %restrict.Omega = Omegai;
    %sh = ones(length(pare),1);
    %if Pbull<0 
    %    restrict.scale = [1./sh];
    %else
    %    restrict.scale = [ones(sizOm,1);1./sh];
    %end
    %pare = pare.*sh;
    paren = pare .*(1 + randn(length(pare),1)*0.001); 
    [pare,fval,exitflag] = fminunc(@(x) cal_quasi_like_BLT(x,z,n,sf,index,Pbull,restrict),paren,options);  
    [llc2n,resc2] = cal_quasi_like_BLT(pare,z,n,sf,index,Pbull,restrict);

    if (llc2n<llc2) % new estimate is better
        pare = paren;
        llc2 = llc2n; 
    end;
    %pare = pare./sh;
end;

[llc2,resc2] = cal_quasi_like_BLT(pare,z,n,sf,index,Pbull,restrict);



