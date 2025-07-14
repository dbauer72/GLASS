function result = random_improve_BLT(result,k,M);
% random_improve_call perturbs the parameter vector in result 
% using a normally distributed disturbance with variance k
% and optimizes the criterion function.
% M random draws are chosen. 
%
% INPUT: result ... result structure
%
% OUTPUT: result ... updated result structure
%
% AUTHOR: dbauer, 4.6.2025


Pbull = -1; 

% extract characteristics
n = result.n;
sf = result.s;
index = result.restrict;
z = result.y; 
pare = result.param; 

[llc2,resc2] = cal_quasi_like_BLT(pare,z,n,sf,index,Pbull);

options = optimoptions('fminunc','display','iter','MaxIterations',1000);
options.MaxFunctionEvaluations = 20000;

for j=1:M

    

    paren = pare + randn(length(pare))*k; 
    [paren,fval,exitflag] = fminunc(@(x) cal_quasi_like_BLT(x,z,n,sf,index,Pbull),paren,options);  
    [llc2n,resc2] = cal_quasi_like_BLT(paren,z,n,sf,index,Pbull);

    disp(sprintf('Running random improve %2d. Old value: %4.2f -> new value: %4.2f.',j,llc2,llc2n));

    if (llc2n<llc2) % new estimate is better
        pare = paren;
        llc2 = llc2n; 
    end;
    %pare = pare./sh;
end;

[llc2,resc2] = cal_quasi_like_BLT(pare,z,n,sf,index,Pbull);
result = compile_results_BLT(pare,n,sf,index,z,Pbull,n);
