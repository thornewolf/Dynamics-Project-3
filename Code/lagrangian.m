clear;close all;clc

% Constants
c.m2 = 1; % kg
c.k = 25; % N/m
c.L0 = 0.5;
c.l2 = 1; % m
c.g = 9.81; % m/s^2

syms l2 m2 L0 g thetat thetadot thetaddot phit phidot phiddot l1t l1dot l1ddot k
syms theta(t) phi(t) l1(t)
pa = [0 0];
pg = [l1*sin(theta)+l2/2*sin(phi) -l1*cos(theta)-l2/2*cos(phi)];
r_ag = pg - pa;
v_g = diff(r_ag,t);
v_g = subs(v_g, [theta, diff(theta,t), phi, diff(phi,t), l1, diff(l1,t)], [thetat, thetadot, phit, phidot, l1t, l1dot]);
v_g = v_g(1);
v_g_squared = v_g(1)*v_g(1)+v_g(2)*v_g(2);
T1 = 0;
T2 = 1/2*m2*v_g_squared + 1/24*m2*l2^2*phidot^2;
T = T1 + T2;
V = 1/2*k*(l1t-L0)^2 - m2*g*(l1t*cos(thetat)+l2/2*cos(phit));
Lagrangian = T - V;

% theta
dLdthetadot = diff(Lagrangian,thetadot);
dLdthetadot_subbed = subs(dLdthetadot, [thetat, thetadot, phit, phidot, l1t, l1dot],...
    [theta, diff(theta,t), phi, diff(phi,t), l1, diff(l1,t)]);
ddLdthetadotdt = diff(dLdthetadot_subbed,t);

dLdtheta = diff(Lagrangian,thetat);

eqn(1) = ddLdthetadotdt - dLdtheta == 0;

% phi
dLdphidot = diff(Lagrangian,phidot);
dLdphidot_subbed = subs(dLdphidot, [thetat, thetadot, phit, phidot, l1t, l1dot],...
   [theta, diff(theta,t), phi, diff(phi,t), l1, diff(l1,t)]);
ddLdphidotdt = diff(dLdphidot_subbed,t);

dLdphi = diff(Lagrangian,phit);

eqn(2) = ddLdphidotdt - dLdphi == 0;

% l1
dLdl1dot = diff(Lagrangian,l1dot);
dLdl1dot_subbed = subs(dLdl1dot, [thetat, thetadot, phit, phidot, l1t, l1dot],...
   [theta, diff(theta,t), phi, diff(phi,t), l1, diff(l1,t)]);
ddLdl1dotdt = diff(dLdl1dot_subbed,t);

dLdl1 = diff(Lagrangian,l1t);

eqn(3) = ddLdl1dotdt - dLdl1 == 0;
% done

eqn(:) = subs(eqn(:), [theta diff(theta,t) diff(theta,t,t) phi diff(phi,t)...
    diff(phi,t,t) l1 diff(l1,t) diff(l1,t,t)],...
    [thetat thetadot thetaddot phit phidot phiddot l1t l1dot l1ddot]);

syms thetadot(t) phidot(t) l1dot(t)
x = solve(eqn,[thetaddot, phiddot, l1ddot]);
thetaEOM = subs(x.thetaddot, [thetat, 'thetadot', phit, 'phidot', l1t, 'l1dot'], [theta, thetadot, phi, phidot, l1, l1dot])
phiEOM = subs(x.phiddot, [thetat, 'thetadot', phit, 'phidot', l1t, 'l1dot'], [theta, thetadot, phi, phidot, l1, l1dot])
l1EOM = subs(x.l1ddot, [thetat, 'thetadot', phit, 'phidot', l1t, 'l1dot'], [theta, thetadot, phi, phidot, l1, l1dot])

EOM = odeFunction([thetadot; thetaEOM; phidot; phiEOM; l1dot; l1EOM], [theta; thetadot; phi; phidot; l1; l1dot],l2, m2, L0, g, k)

[T,S] = ode45(@(t,s)EOM(t,s,c.l2, c.m2, c.L0, c.g, c.k), linspace(0,10,1001), [0, 0, 0, 0, c.L0, 0]);
figure(1)
plot(T,S(:,5))

[T,S] = ode45(@(t,s)EOM(t,s,c.l2, c.m2, c.L0, c.g, c.k), linspace(0,10,1001), [pi/18, 0, pi/9, 0, c.L0, 0]);
figure(2)
plot(T,S(:,5))

[T,S] = ode45(@(t,s)EOM(t,s,c.l2, c.m2, c.L0, c.g, c.k), linspace(0,10,1001), [pi/6, 0, pi/3, 0, c.L0, 0]);
figure(3)
plot(T,S(:,5))
