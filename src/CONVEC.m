function [u1,v1,Fu_nm,Fv_nm]=CONVEC(u0,v0,Re,NX,NY,dx,dy,dt,Fu_nm,Fv_nm,iter)

% This function solves the non-linear parts of the 2D Navier-Stokes
% equation. The second order central difference discretization is used for
% spacial derivatives.

% Calculate Fu and Fv
for j=2:NX-1
    jm = j-1;
    jp = j+1;
    for k=2:NY-1
        km = k-1;
        kp = k+1;
        Fu(k,j) = -u0(k,j)*(u0(k,jp)-u0(k,jm))/(2*dx)-v0(k,j)*(u0(kp,j)-u0(km,j))/(2*dy)+(3/Re);        
        Fv(k,j) = -u0(k,j)*(v0(k,jp)-v0(k,jm))/(2*dx)-v0(k,j)*(v0(kp,j)-v0(km,j))/(2*dy);
    end
end

% Set Fu and Fv for FIRST iteration only
if iter == 2
    Fu_nm = Fu;
    Fv_nm = Fv;
end

% Calculate u* and v*
u1(2:NY-1,2:NX-1) = dt*((3/2)*Fu(2:NY-1,2:NX-1) - (1/2)*Fu_nm(2:NY-1,2:NX-1)) + u0(2:NY-1,2:NX-1);
v1(2:NY-1,2:NX-1) = dt*((3/2)*Fv(2:NY-1,2:NX-1) - (1/2)*Fv_nm(2:NY-1,2:NX-1)) + v0(2:NY-1,2:NX-1);

% Update Fu and Fv for next iteration
Fu_nm = Fu;
Fv_nm = Fv;
end


