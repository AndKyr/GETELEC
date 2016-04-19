function [J,G]=J_sph_approx(F,W,R)
dW=0.001;


G=G_approx(F,W,R);

if exp(-G)>0.5
    J=1;
    return;
end
G2=G_approx(F,W+dW,R);

d=abs(G2-G)/dW;
J=1.62e-4*(d.^-2).*exp(-G);
end

function G=G_approx(F,W,R)

Vfun=@(x) W-F*R*x./(x+R)-0.36./(x+0.5*x.^2/R);
Vfun2= @(x) sqrt(W-F*R*x./(x+R)-0.36./(x+0.5*x.^2/R));

FNpoly=[-F W -0.36];
FNzeros=roots(FNpoly);
x0=min(FNzeros);
x1=max(FNzeros)+0.05;
if abs(x1-x0)<0.1
    G=0;
    return;
end

x0=fzero(Vfun,x0);
x1=fzero(Vfun,x1);

G=10.246*quadgk(Vfun2,x0+1e-6,x1-1e-6);
end