function [td]=VISCsolver(te,jmax,dx,dt,Re)
% This function solves the 1D diffustion problem. The functino is a
% modified version of "difim_mat.m" used for AME535a homework assignment.

% Parameters
alph=1/Re;
s=alph*dt/(dx^2);
gam=0;
jmap=jmax-1;
bet = 0.5 + gam;
emx(1)=0;
emx(2)=1-(2*emx(1));
emx(3)=emx(1);
ad=(emx(1)*(1+gam))-(bet*s);
bd=(emx(2)*(1+gam))+(2*bet*s);
cd=ad;
elx(1)=1;
elx(2)=-2;
elx(3)=1;
jmaf=jmax-2;
n=1;

% Use exact solution as initial condition

tol = te;
tn  = te;

% set up the tridiagonal system of equations
for j=2:jmap
    jm=j-1;
    if(n<=1)
        A(1,jm)=0;
        A(2,jm)=ad;
        A(3,jm)=bd;
        A(4,jm)=cd;
        A(5,jm)=0;
    end
    
    d(jm)=0;
    for k=1:3
        kj=j-2+k;
        d(jm)=d(jm)+(emx(k)*(((1+(2*gam))*tn(kj))-(gam*tol(kj))));
        d(jm)=d(jm)+(s*elx(k)*(1-bet)*tn(kj));
    end
      
end

d(1)=d(1)-(A(2,1)*tn(1));
d(jmaf)=d(jmaf)-(A(4,jmaf)*tn(jmax));

% solve banded system of equations

if(n==1)
    A=banfac(A,jmaf);
end

dum=bansol(d,A,jmaf);

for j=2:jmap
    tol(j)=tn(j);
    tn(j)=dum(j-1);
end

for j=2:jmap
    td(j)=tn(j);
end

td(jmax) = te(jmax);
end
