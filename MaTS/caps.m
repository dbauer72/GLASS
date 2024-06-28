function x = caps(x,M);
% caps off the values in x at plus minus M.
%

if nargin<2
    M= 10^6;
end

x(x>M)=M+ (exp(x(x>M)/M-1)-1)/(exp(x(x>M)/M-1)+1);
x(x<-M)=-M- (exp(-x(x>M)/M-1)-1)/(exp(-x(x>M)/M-1)+1);
