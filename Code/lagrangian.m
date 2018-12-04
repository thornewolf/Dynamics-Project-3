clear;close all;clc

% Constants
c.m2 = 1; % (kg) bar mass
c.k = 25; % (N/m) spring constant
c.L0 = 0.5; % (m) spring unstretched length
c.l2 = 1; % (m) length of bar
c.g = 9.81; % (m/s^2) gravity

syms l2 m2 L0 g thetat thetadot thetaddot phit phidot phiddot l1t l1dot l1ddot k
syms theta(t) phi(t) l1(t)

% Defining origin as point A
pa = [0 0];
% Defining center of mass of bar as point G
pg = [l1*sin(theta)+l2/2*sin(phi) -l1*cos(theta)-l2/2*cos(phi)];
% Vector to G from A
r_ag = pg - pa;
% Time derivative of vector is velocity
v_g = diff(r_ag,t);
v_g = subs(v_g, [theta, diff(theta,t), phi, diff(phi,t), l1, diff(l1,t)],...
                [thetat, thetadot, phit, phidot, l1t, l1dot]);
v_g = v_g(1);
v_g_squared = v_g(1)*v_g(1)+v_g(2)*v_g(2);

% Kinetic Energy of Spring
T1 = 0;
% Kinetic Energy of Bar
T2 = 1/2*m2*v_g_squared + 1/24*m2*l2^2*phidot^2;
% Total Kinetic Energy
T = T1 + T2;

% Potential Energy of Spring
V1 = 0;
% Potential Energy of Bar
V2 = 1/2*k*(l1t-L0)^2 - m2*g*(l1t*cos(thetat)+l2/2*cos(phit));
% Total Potential Energy
V = V1 + V2;

Lagrangian = T - V;

% Partial of Lagrange Eq. w.r.t. thetadot
dLdthetadot = diff(Lagrangian,thetadot);
dLdthetadot_subbed = subs(dLdthetadot, [thetat, thetadot, phit, phidot, l1t, l1dot],...
    [theta, diff(theta,t), phi, diff(phi,t), l1, diff(l1,t)]);
% Time derivative of Partial of Lagrange Eq. w.r.t. thetadot
ddLdthetadotdt = diff(dLdthetadot_subbed,t);
% Partial of Lagrange Eq. w.r.t. theta
dLdtheta = diff(Lagrangian,thetat);

% Theta EOM 
eqn(1) = ddLdthetadotdt - dLdtheta == 0;

% Partial of Lagrange Eq. w.r.t. phidot
dLdphidot = diff(Lagrangian,phidot);
dLdphidot_subbed = subs(dLdphidot, [thetat, thetadot, phit, phidot, l1t, l1dot],...
   [theta, diff(theta,t), phi, diff(phi,t), l1, diff(l1,t)]);
% Time derivative of Partial of Lagrange Eq. w.r.t. phidot
ddLdphidotdt = diff(dLdphidot_subbed,t);
% Partial of Lagrange Eq. w.r.t. phi
dLdphi = diff(Lagrangian,phit);

% Phi EOM
eqn(2) = ddLdphidotdt - dLdphi == 0;

% Partial of Lagrange Eq. w.r.t. l (spring deflection)
dLdl1dot = diff(Lagrangian,l1dot);
dLdl1dot_subbed = subs(dLdl1dot, [thetat, thetadot, phit, phidot, l1t, l1dot],...
   [theta, diff(theta,t), phi, diff(phi,t), l1, diff(l1,t)]);
% Time derivative of Partial of Lagrange Eq. w.r.t. l
ddLdl1dotdt = diff(dLdl1dot_subbed,t);
% Partial of Lagrange Eq. w.r.t. l
dLdl1 = diff(Lagrangian,l1t);

% l (spring) EOM 
eqn(3) = ddLdl1dotdt - dLdl1 == 0;

eqn(:) = subs(eqn(:), [theta diff(theta,t) diff(theta,t,t) phi diff(phi,t)...
    diff(phi,t,t) l1 diff(l1,t) diff(l1,t,t)],...
    [thetat thetadot thetaddot phit phidot phiddot l1t l1dot l1ddot]);

syms thetadot(t) phidot(t) l1dot(t)

x = solve(eqn,[thetaddot, phiddot, l1ddot]);
thetaEOM = subs(x.thetaddot, [thetat, 'thetadot', phit, 'phidot', l1t, 'l1dot'], ...
                             [theta, thetadot, phi, phidot, l1, l1dot]);
                         
phiEOM = subs(x.phiddot, [thetat, 'thetadot', phit, 'phidot', l1t, 'l1dot'],...
                         [theta, thetadot, phi, phidot, l1, l1dot]);
                     
l1EOM = subs(x.l1ddot, [thetat, 'thetadot', phit, 'phidot', l1t, 'l1dot'], ...
                       [theta, thetadot, phi, phidot, l1, l1dot]);

EOM = odeFunction([thetadot; thetaEOM; phidot; phiEOM; l1dot; l1EOM],...
                  [theta; thetadot; phi; phidot; l1; l1dot],l2, m2, L0, g, k);

              
[T,S] = ode45(@(t,s)EOM(t,s,c.l2, c.m2, c.L0, c.g, c.k), linspace(0,10,1001),...
    [0, 0, 0, 0, c.L0, 0]);
figure(1)
hold on
grid on
xlabel('Time, sec')
ylabel('Length, m')
plot(T,S(:,5),'DisplayName', 'Spring Length, r')
legend('show','Location','northoutside')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 6 2.5]);
print('0-0','-dpdf');
[T,S] = ode45(@(t,s)EOM(t,s,c.l2, c.m2, c.L0, c.g, c.k), linspace(0,10,1001),...
    [pi/18, 0, pi/9, 0, c.L0, 0]);
figure(2)
hold on
grid on
xlabel('Time, sec')
ylabel('Length, m')
plot(T,S(:,5),'DisplayName', 'Spring Length, r')
legend('show','Location','northoutside')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 6 2.5]);
print('18-9:Spring','-dpdf');
figure(3)
hold on
grid on
xlabel('Time, sec')
ylabel('Angular Position, rad')
plot(T,S(:,1),'DisplayName', 'Spring Angular Deflection, \theta')
plot(T,S(:,3),'DisplayName', 'Bar Angular Deflection, \phi')
legend('show','Location','northoutside')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 6 2.5]);
print('18-9:Angles','-dpdf');

[T,S] = ode45(@(t,s)EOM(t,s,c.l2, c.m2, c.L0, c.g, c.k), linspace(0,10,1001),...
    [pi/6, 0, pi/3, 0, c.L0, 0]);
figure(4)
hold on
grid on
xlabel('Time, sec')
ylabel('Length, m')
plot(T,S(:,5),'DisplayName', 'Spring Length, r')
legend('show','Location','northoutside')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 6 2.5]);
print('6-3:Spring','-dpdf');
figure(5)
hold on
grid on
xlabel('Time, sec')
ylabel('Angular Position, rad')
plot(T,S(:,1),'DisplayName', 'Spring Angular Deflection, \theta')
plot(T,S(:,3),'DisplayName', 'Bar Angular Deflection, \phi')
legend('show','Location','northoutside')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 6 2.5]);
print('6-3:Angles','-dpdf');
