function s=cprint(n,x)
%
if n<=8
   error('n must be larger than 8!');
end

ax=abs(x);
l=(10^(n-2)-0.5)*10.^((-n+2:1)');

if ax>=(10^(n-1)-0.5)
   fmt=sprintf('%%%i.%iE',n,n-8);
elseif ax>=(10^(n-2)-0.5)
   fmt=sprintf('%%%i.0f',n);
elseif ax>=l(1)
   d=sum(ax>=l);
   fmt=sprintf('%%%i.%if',n,n-d-2);
elseif (ax<(2^(n-7)-0.5)^(-6))&(ax>0)
   fmt=sprintf('%%%i.%iE',n,n-8);
else
   fmt=sprintf('%%%i.%if',n,n-3);
end
s=sprintf(fmt,x);
