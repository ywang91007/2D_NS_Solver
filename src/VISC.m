function [u4,v4]=VISC(u2,v2,Re,NX,NY,dx,dy,dt)
% This function solves the viscous terms in the 2D Navier-Stokes equation
% Sub-function "VISCsolver.m" is used to calculate the viscous term


% preallocation
u3=u2;
v3=v2;

% Viscous term in x direction
for k=1:NY
    u3(k,:)=VISCsolver(u2(k,:),NX,dx,dt/2,Re);
    v3(k,:)=VISCsolver(v2(k,:),NX,dx,dt/2,Re);
end


% preallocation
u4=u3;
v4=v3;

% Viscous term in y direction
for j=1:NX
    u4(:,j)=VISCsolver(u3(:,j),NY,dy,dt/2,Re);
    v4(:,j)=VISCsolver(v3(:,j),NY,dy,dt/2,Re);
end

% Boundary conditions

% Bottom (No slip)
u4(1,:) = 0;
v4(1,:) = 0;

% Top (No slip)
u4(NY,:) = 0;
v4(NY,:) = 0;

% Inlet (uniform)
u4(:,1) = 1;
v4(:,1) = 0;

% Outlet (Neumann condition)
u4(:,NX) = u4(:,NX-1);
v4(:,NX) = 0;

end
