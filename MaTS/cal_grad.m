function [crit,grad,Hess,Jgrad,JHess,ve,psi] = cal_grad(Y,theta,L,p);
% evaluates gradient and Hessian numerically to compare with analytic
% expression. 


[crit,grad,Hess,ve,psi] = cal_crit_MaTS(Y,theta,L,p);

Jgrad = grad*0;
JHess = Hess*0;

% gradient: 
eps = 0.00001;
for j=1:length(theta)
    thetan = theta;
    thetan(j) = thetan(j)+ eps;
    crn = cal_crit_MaTS(Y,thetan,L,p);
    Jgrad(j) = (crn-crit)/eps;
end

% Hessian: four point approximation 
for jr=1:length(theta)
    % two points in direction jr. 
    thetar = theta;
    thetar(jr)=thetar(jr)+eps; 
    thetarm = theta;
    thetarm(jr)=thetarm(jr)-eps; 
    for jc = 1:jr
        % four points in direction jc.
        thetapp = thetar;
        thetapp(jc)=thetapp(jc)+eps; 

        thetapm = thetar;
        thetapm(jc)=thetapm(jc)-eps; 

        thetamp = thetarm;
        thetamp(jc)=thetamp(jc)+eps; 

        thetamm = thetarm;
        thetamm(jc)=thetamm(jc)-eps; 

        % P1: pp
        crpp = cal_crit_MaTS(Y,thetapp,L,p);
        % P2
        crpm = cal_crit_MaTS(Y,thetapm,L,p);

        % P3
        crmp = cal_crit_MaTS(Y,thetamp,L,p);

        %P4
        crmm = cal_crit_MaTS(Y,thetamm,L,p);

        % calculate
        JHess(jr,jc) = (crpp - crmp - crpm + crmm)/(4*eps*eps);
    end
end
dJH = diag(JHess);
JHess = JHess + JHess'-diag(dJH);
