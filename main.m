% Author: Yi Wang
% Date: 2024.12.12

% This is the main script that calculate the 2D Navier-Stokes equations for
% a steady, viscous entarance flow into a 2D plane channel.
% Sub-fuction "Jacobian.m" is used to calculate the Jacobian matrix used
% for pressure calculation
% Sub-function "CONVEC.m" is used to calculate the non-linear terms
% Sub-function "PRESS.m" is used to calulate the pressure term
% Sub-function "VISC.m" is ued to calulate the viscous terms

clear all; close all; clc;

% User defined parameters
Re = 10;
NX = 51;
NY = 51;
dt = 0.01;
tmax = 20;
tolerance = 0.001;
xmax = 10;
ymax = 2;


% Variables
x = linspace(0,xmax,NX);
y = linspace(-1,1,NY);
dx = xmax/(NX-1);
dy = ymax/(NY-1);
n = NX*NY;
Fu_nm = 0;
Fv_nm = 0;
time = 0;
iter = 1;
diff = ones(1,tmax/dt);
RHS = zeros(NY,NX);
Q = zeros(1,NX);

% Initial conditions (uniform u0, zero v0)
u0 = ones(NY,NX);
v0 = zeros(NY,NX);
u0(1,:) = 0;
u0(NY,:) = 0;
u1 = u0;
v1 = v0;

% Exact solution at outlet
u_exact = (3/2)*(1-y.^2);
u_exact = u_exact.';

% Jacobian matrix for pressure calculation
[A,jpvt] = Jacobian(NX,NY,dx,dy);

% Iterations
while diff(iter) > tolerance && time < tmax
    iter = iter+1;
    % Nonlinear step
    [u1,v1,Fu_nm,Fv_nm] = CONVEC(u0,v0,Re,NX,NY,dx,dy,dt,Fu_nm,Fv_nm,iter);    
    % Pressure step
    [u2,v2,P] = PRESS(u1,v1,NX,NY,dx,dy,dt,A,jpvt,n,RHS);
    % Viscous step
    [u4,v4] = VISC(u2,v2,Re,NX,NY,dx,dy,dt);
    
    diff(iter) = sum(sum(abs(u4-u0)+abs(v4-v0)));
    u_rms(iter) = sqrt(sum((u_exact-u4(:,NX)).^2/(NY-2)));
    Q_rms(iter) = sqrt(sum((2-sum(u4(:,NX))*dy).^2/NX));
    u0 = u4;
    v0 = v4;
    time = time+dt;
end

iteration = 1:iter;
% Final results
u_n = u4;
v_n = v4;

% show rms error and iterations
disp(['u_rms = ',num2str(u_rms(end))])
disp(['Q_rms = ',num2str(Q_rms(end))])
disp(['iterations = ',num2str(iter)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j_vals = [2, (NX+1)/2, NX-1, NX];

% Plotting u-profiles
figure();
plot(u_n(:,j_vals(1)), y, 'b-', ...            % profile at j=2
     u_n(:,j_vals(2)), y, 'r-', ...            % profile at j=(nx+1)/2
     u_n(:,j_vals(3)), y, 'g-', ...            % profile at j=nx-1
     u_n(:,j_vals(4)), y, 'k--')               % profile at j=nx
hold on;
plot(u_exact, y, 'm-', 'LineWidth', 1.5);
xlabel('u');          
ylabel('y');          
title('Vertical Profiles of u at Selected x-locations');
legend('j=2', ['j=' num2str((NX+1)/2)], ['j=' num2str(NX-1)], ['j=' num2str(NX)], 'Exact outflow', ...
       'Location', 'best');
grid on;

% Plotting v-profiles
figure();
plot(v_n(:,j_vals(1)), y, 'b-', ...
     v_n(:,j_vals(2)), y, 'r-', ...
     v_n(:,j_vals(3)), y, 'g-', ...
     v_n(:,j_vals(4)), y, 'k--')

xlabel('v');
ylabel('y');
title('Vertical Profiles of v at Selected x-locations');
legend(['j=2'], ['j=' num2str((NX+1)/2)], ['j=' num2str(NX-1)], ['j=' num2str(NX)], ...
       'Location', 'best');
grid on;

% Plot u at centerline
kc = (NY+1)/2;
u_center = u_n(kc, :);  
v_center = v_n(kc, :);  
figure();
plot(x, u_center, 'b-', 'LineWidth', 1.5);
hold on;
yline(1.5, 'r--', 'Asymptotic u=1.5');
xlabel('x');
ylabel('u');
title('Velocity u along the horizontal centerline');
legend('u(x) at y_{center}', 'Asymptotic Value', 'Location', 'best');
grid on;

% Plot v at centerline
figure();
plot(x, v_center, 'g-', 'LineWidth', 1.5);
xlabel('x');
ylabel('v');
title('Velocity v along the horizontal centerline');
grid on;

% Plot the flow rate Q(x)
for j = 1:NX
    Q(j) = sum(u_n(:,j))*dy;  
end
Q_inflow = sum(u_n(:,1))*dy;
figure();
plot(x, Q, 'b-', 'LineWidth', 1.5);
hold on;
yline(Q_inflow, 'r--', 'Inflow Rate', 'LineWidth', 1.5);
xlabel('x');
ylabel('Q(x)');
title('Flow Rate Q(x) along the Channel');
legend('Q(x)', 'Inflow Rate', 'Location', 'best');
grid on;

figure()
plot(iteration, u_rms)