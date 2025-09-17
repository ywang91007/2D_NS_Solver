function [A,jpvt]=Jacobian(NX,NY,dx,dy)
% This function calculate the Jacobian matrix used for pressure calculation
% Sub-function "fact.m" is used to apply LU decompostion to the Jacobian
% matrix.

% Define coefficients
n=NX*NY;
dj=linspace(1,NX,NX);
dk=linspace(1,NY,NY);
E1=1/dx^2;
E2=1/dy^2;
E3=-2/dx^2-2/dy^2;
A=zeros(n);
jpvt = zeros(1,n);

% Dummy variables
a=zeros(n);
b=a;
p=1;
q=1;

% Jacobian calculation
while p<=n
    a(p:p+NY-1,:)=dj(q);
    b(p:p+NY-1,1)=dk;
    p=p+NY;
    q=q+1;
end
b(:,2:end) = repmat(b(:,1),1,size(b,2)-1);
c=a.';
d=b.';

for j=1:n
    for k=1:n
        if (j==k)
            A(k,j)=E3;
        end
        if (a(k,j)==c(k,j)) && (abs(b(k,j)-d(k,j))==1)
            A(k,j)=E2;
        end
        if (b(k,j)==d(k,j)) && (abs(a(k,j)-c(k,j))==1)
            A(k,j)=E1;
        end
    end
end

% LU decomposition
[A,jpvt]=fact(n,A,jpvt);

end