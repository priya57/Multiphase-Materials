% Newmarks Method
clear;
clc;

%% Initialization
nEle = 1;
nNod = nEle+1;
nDof = nNod*2;
E = 1000;
mu = 0;
rho = 1000;
A = 1;
L = 1.0/nEle;
kDarcy = 0;

%% Initializing the matrices
stiff = zeros(nDof,nDof);
damp = zeros(nDof,nDof);
mass = zeros(nDof,nDof);

k1 = zeros(nDof,nDof);
k2 = zeros(nDof,nDof);
k3 = zeros(nDof,nDof);

%% Calculation of the mass, damping and stifness matrices
for i = 1:2:nDof-2
    mass(i,i)       = mass(i,i) + rho*A*L/3;
    mass(i,i+2)     = rho*A*L/6;
    mass(i+2,i)     = rho*A*L/6;
    mass(i+2,i+2)   = rho*A*L/3;
    
    damp(i,i)       = damp(i,i) + mu*A/L;
    damp(i,i+2)     = -mu*A/L;
    damp(i+2,i)     = -mu*A/L;
    damp(i+2,i+2)   = mu*A/L;
    
    stiff(i,i)       = stiff(i,i) + E*A/L;
    stiff(i,i+2)     = -E*A/L;
    stiff(i+2,i)     = -E*A/L;
    stiff(i+2,i+2)   = E*A/L;
    
    k1(i+1,i+1) = k1(i,i) + kDarcy*A/L;
    k1(i+1,i+3) = -kDarcy*A/L;
    k1(i+3,i+1) = -kDarcy*A/L;
    k1(i+3,i+3) = kDarcy*A/L;
    
    k2(i+1,i)   = k2(i,i) - 0.5*A;
    k2(i+1,i+2) = 0.5*A;
    k2(i+3,i)   = -0.5*A;
    k2(i+3,i+2) = 0.5*A;
    
    k3(i,i+1)   = k3(i,i) + 0.5*A;
    k3(i,i+3)   = 0.5*A;
    k3(i+2,i+1) = -0.5*A;
    k3(i+2,i+3) = -0.5*A;
end

%% Initializing time steps
T = 100;            % Total time
dt = 0.1;           % Time increment
t = 0:dt:T;         % Time step array
n = numel(t);       % Number of time steps

%% Initializing force matrix
f = 0;             % Maximum force
omega = 1;
F = zeros(nDof,n);
F(nDof,:) = f*cos(omega.*t);

% Solving the matrices by Newmark's method
%% Initializing the displacement, velocity and acceleration matrices
dx = zeros(nDof,n);
vdisp = zeros(nDof,n);
vel = zeros(nDof,n);
acc = zeros(nDof,n);

%% Defining the intial values of displacement, velocity and acceleration
vdisp((nDof-1):nDof,1) = 0.1;    % initial displacement
vel(nDof,1) = 0;                % initial velocity
acc(nDof,1) = 0;                % initial acceleration

%% Global stiffness matrix
K = stiff + k1+ k2 + k3;

%% LHS Matrix
LHS = mass*4/dt^2 + damp*2/dt + K;
LHSred = LHS(3:nDof,3:nDof);

%% RHS Matrix
for i = 2:n
    RHS = F(:,i) + mass*acc(:,i-1) + (4/dt*mass + damp)*vel(:,i-1) - K*vdisp(:,i-1);
    RHSred = RHS(3:nDof);
   
    dxred = LHSred\RHSred;
    dx = [0; 0; dxred];
    
    vdisp(:,i) = vdisp(:,i-1) + dx;
    vel(:,i) = 2/dt*dx - vel(:,i-1);
    acc(:,i) = 4/dt^2*dx - 4/dt*vel(:,i-1) - acc(:,i-1);   
end

 disp= zeros(nNod,n);
for i = 1:nNod
    disp(i,:) = vdisp(2*i-1,:);
end

%% Plotting
figure
plot(t,disp)
xlabel('Time [s]')
ylabel('Displacement [m]')

figure
plot(t,vdisp)
xlabel('Time [s]')
ylabel('Vibrational Displacement [m]')

figure
plot(t,vel)
xlabel('Time [s]')
ylabel('Velocity [m/s]')

figure
plot(t,acc)
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')