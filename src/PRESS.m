function [u2,v2,P] = PRESS(u1,v1,NX,NY,dx,dy,dt,A,jpvt,n,RHS)
% This function solves the pressure term in the 2D Navier-Stokes equation.
% Sub-function "solve.m" is used to calculate the pressure.

% Calculate RHS
for j=2:NX-1
    jm = j-1;
    for k=2:NY-1
        km = k-1;
        RHS(k,j) = (1/dt)*((u1(k,j)-u1(k,jm))/dx+(v1(k,j)-v1(km,j))/dy);
    end
end

% Calculate pressure
P = solve(n-1,A,jpvt,RHS);   

% Calculate pressures at four corners
P(1,1) = (P(1,2)+P(2,1))/2;
P(1,NX) =(P(1,NX-1)+P(2,NX))/2;
P(NY,1) = (P(NY,2)+P(NY-1,1))/2;
P(NY,NX) = (P(NY,NX-1)+P(NY-1,NX))/2;

% Calculate u** and v**
for j=2:NX-1
    jp = j+1;
    jm = j-1;
    for k=2:NY-1
        kp = k+1;
        km = k-1;
        u2(k,j) = u1(k,j)-dt*((P(k,jp)-P(k,jm))/(2*dx));
        v2(k,j) = v1(k,j)-dt*((P(kp,j)-P(km,j))/(2*dy));
    end
end

% Boundary conditions

% Bottom (No slip)
u2(1,:) = 0;
v2(1,:) = 0;

% Top (No slip)
u2(NY,:) = 0;
v2(NY,:) = 0;

% Left (uniform)
u2(:,1) = 1;
v2(:,1) = 0;

% Right (Neumann condition)
u2(:,NX) = u2(:,NX-1);
v2(:,NX) = 0;

end




