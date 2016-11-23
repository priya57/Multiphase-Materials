%% Newmarks Method
clear;
clc;

%% Initialization
nEle = 1;
nNod = nEle+1;
E = 1000;
mu = 10;
rho_0 = 1000;
A0= 1;
L0= 1.0/nEle;

%% Initializing the matrices
stiff = zeros(nNod,nNod);
damp = zeros(nNod,nNod);
mass = zeros(nNod,nNod);

%% Calculation of the mass, damping and stifness matrices
for i = 1:nNod
    for j = 1:nNod
        if i == j && (i == 1 || i == nNod)
            stiff(i,j) = E*A0/L0;
            damp(i,j) = mu*A0/L0;
            mass(i,j) = rho_0*A0*L0/3;
        elseif i == j && (i ~= 1 || i ~= nNod)
            stiff(i,j) = 2*E*A0/L0;
            damp(i,j) = 2*mu*A0/L0;
            mass(i,j) = 2*rho_0*A0*L0/3;
        
        elseif i == j+1 || j == i+1
            stiff(i,j) = -E*A0/L0;
            damp(i,j) = -mu*A0/L0;
            mass(i,j) = rho_0*A0*L0/6;
        end
    end
end

%Defining the boundary conditions for the first node to be fixed
%stiff(1,:) = 0;
%stiff(:,1) = 0;
%damp(1,:) = 0;             % misstake
%damp(:,1) = 0;
%mass(1,:) = 0;
%mass(:,1) = 0;

%% Initializing time steps
T = 100;             % Total time
dt = 0.1;           % Time increment
t = 0:dt:T;         % Time step array
n = numel(t);       % Number of time steps

%% Initializing force matrix
f = 0;             % Maximum force
omega = 1;
F = zeros(nNod,n);
% Effect of force on each node
% (boundary condition for the first node is already included)
%F(nNod-1,:) = (f/2)*cos(omega.*t);
F(nNod,:) = f*cos(omega.*t);
%F(2,:) = f*cos(omega.*t);

%% Solving the matrices by Newmark's method

% Initializing the displacement, velocity and acceleration matrices
disp = zeros(nNod,n);
vel = zeros(nNod,n);
acc = zeros(nNod,n);

%dx = zeros(nNod,n);
%% Calculating the intial values of displacement, velocity and acceleration
disp(nNod,1) = 0.1;
vel(nNod,1) = 0;
acc(nNod,1) = 0;

%dx(:,1) = (mass.*d2t + damp.*d1t + stiff)\F(:,1);
%vel(:,1)= d1t*dx(:,1);
%acc(:,1)= d2t*dx(:,1);
%dx(1,:) = 0;
%% LHS Matrix

LHS = mass*4/dt^2 + damp*2/dt + stiff;
LHSred=LHS;
LHSred = LHS(2:nNod,2:nNod);

%% RHS Matrix

for i = 2:n
    RHS = F(:,i) + mass*acc(:,i-1) + (4/dt*mass + damp)*vel(:,i-1) - stiff*disp(:,i-1);
    RHSred = RHS(2:nNod);

    %dxred = dx(nNod-1,nNod);
    dxred = LHSred\RHSred;
    dx = [0; dxred];
    
    %dx(:,i) = [0; dxred(:,i)];
    disp(:,i) = disp(:,i-1) + dx;
    vel(:,i) = 2/dt*dx - vel(:,i-1);
    acc(:,i) = 4/dt^2*dx - 4/dt*vel(:,i-1) - acc(:,i-1);   
end

%% Plotting
figure
plot(t,disp)
xlabel('Time [s]')
ylabel('Displacement [m]')
%legend('Node1','Node2','Node3')

figure
plot(t,vel)
xlabel('Time [s]')
ylabel('Velocity [m/s]')
%legend('Node1','Node2','Node3')

figure
plot(t,acc)
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')
%legend('Node1','Node2','Node3')