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
V1 = 1/2*k*(l1t-L0)^2;
% Potential Energy of Bar
V2 = - m2*g*(l1t*cos(thetat)+l2/2*cos(phit));
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


options = odeset('AbsTol', 1e-9, 'RelTol', 1e-6);

%-------------------------------------------------------------------------%
% 0 and 0
%-------------------------------------------------------------------------%
[T1,S1] = ode45(@(t,s)EOM(t,s,c.l2, c.m2, c.L0, c.g, c.k), linspace(0,10,1001),...
    [0, 0, 0, 0, c.L0, 0]);
[T2,S2] = ode45(@(t,s)EOM(t,s,c.l2, c.m2, c.L0, c.g, c.k), linspace(0,10,1001),...
    [0, 0, 0, 0, c.L0, 0],options);
[T3,S3] = ode45(@(t,s)EOM(t,s,c.l2, c.m2, c.L0, c.g, c.k), linspace(0,10,1001),...
    [pi/18, 0, pi/9, 0, c.L0, 0]);

%-------------------------------------------------------------------------%
[T4,S4] = ode45(@(t,s)EOM(t,s,c.l2, c.m2, c.L0, c.g, c.k), linspace(0,10,1001),...
    [pi/18, 0, pi/9, 0, c.L0, 0],options);

%-------------------------------------------------------------------------%
%  pi/6 and pi/3
%-------------------------------------------------------------------------%

[T5,S5] = ode45(@(t,s)EOM(t,s,c.l2, c.m2, c.L0, c.g, c.k), linspace(0,10,1001),...
    [pi/6, 0, pi/3, 0, c.L0, 0]);

%-------------------------------------------------------------------------%
[T6,S6] = ode45(@(t,s)EOM(t,s,c.l2, c.m2, c.L0, c.g, c.k), linspace(0,10,1001),...
    [pi/6, 0, pi/3, 0, c.L0, 0],options);
syms thetadot phidot l1dot
v_g = diff(r_ag,t)
v_g = subs(v_g, [theta, diff(theta,t), phi, diff(phi,t), l1, diff(l1,t)],...
                [thetat, thetadot, phit, phidot, l1t, l1dot]);
v_g = v_g(1);
SM = [S1, S2, S3, S4, S5, S6];
T = T1;
S = S1;
energy_vector1 = zeros(length(S));
for i = 1:length(S)
    if mod(i,100)==0
        i
    end
    v_g_subbed = eval(subs(v_g, [thetat, thetadot, phit, phidot, l1t, l1dot],...
        [S(i,1), S(i,2), S(i,3), S(i,4), S(i,5), S(i,6)]));
    v_g_squared = v_g_subbed(1)*v_g_subbed(1)+v_g_subbed(2)*v_g_subbed(2);
    integration_Kinetic_Energy = 1/2*c.m2*v_g_squared + 1/24*c.m2*c.l2^2*S(i,4)^2;
    integration_Potential_Energy = 1/2*c.k*(S(i,5)-c.L0)^2 - c.m2*c.g*(S(i,5)*cos(S(i,1))+c.l2/2*cos(S(i,3)));
    total_Energy = integration_Kinetic_Energy + integration_Potential_Energy;
    energy_vector1(i) = total_Energy;
end
energy_vector1 = energy_vector1(:,1);

T = T2;
S = S2;
energy_vector2 = zeros(length(S));
for i = 1:length(S)
    if mod(i,100)==0
        i
    end
    v_g_subbed = eval(subs(v_g, [thetat, thetadot, phit, phidot, l1t, l1dot],...
        [S(i,1), S(i,2), S(i,3), S(i,4), S(i,5), S(i,6)]));
    v_g_squared = v_g_subbed(1)*v_g_subbed(1)+v_g_subbed(2)*v_g_subbed(2);
    integration_Kinetic_Energy = 1/2*c.m2*v_g_squared + 1/24*c.m2*c.l2^2*S(i,4)^2;
    integration_Potential_Energy = 1/2*c.k*(S(i,5)-c.L0)^2 - c.m2*c.g*(S(i,5)*cos(S(i,1))+c.l2/2*cos(S(i,3)));
    total_Energy = integration_Kinetic_Energy + integration_Potential_Energy;
    energy_vector2(i) = total_Energy;
end
energy_vector2 = energy_vector2(:,1);


figure(1)
hold on

grid on

(max(energy_vector1)-max(energy_vector2))/max(energy_vector2)*100

[min(energy_vector1), max(energy_vector1)]
ylim([min(energy_vector1), max(energy_vector1)])
plot(T,energy_vector1,'DisplayName', 'Default Tol')
plot(T,energy_vector2,'DisplayName', 'Increased Tol')
legend('show','Location','northoutside')
xlabel('Time, sec')
ylabel('Total Energy, J')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 6 2.5]);
print('-f1','Energy1','-dpdf');

hold off

T = T3;
S = S3;
energy_vector3 = zeros(length(S));
for i = 1:length(S)
    if mod(i,100)==0
        i
    end
    v_g_subbed = eval(subs(v_g, [thetat, thetadot, phit, phidot, l1t, l1dot, l2],...
        [S(i,1), S(i,2), S(i,3), S(i,4), S(i,5), S(i,6), c.l2]));
    v_g_squared = v_g_subbed(1)*v_g_subbed(1)+v_g_subbed(2)*v_g_subbed(2);
    integration_Kinetic_Energy = 1/2*c.m2*v_g_squared + 1/24*c.m2*c.l2^2*S(i,4)^2;
    integration_Potential_Energy = 1/2*c.k*(S(i,5)-c.L0)^2 - c.m2*c.g*(S(i,5)*cos(S(i,1))+c.l2/2*cos(S(i,3)));
    total_Energy = integration_Kinetic_Energy + integration_Potential_Energy;
    energy_vector3(i) = total_Energy;
end
energy_vector3 = energy_vector3(:,1);

T = T4;
S = S4;
energy_vector4 = zeros(length(S));
for i = 1:length(S)
    if mod(i,100)==0
        i
    end
    v_g_subbed = eval(subs(v_g, [thetat, thetadot, phit, phidot, l1t, l1dot, l2],...
        [S(i,1), S(i,2), S(i,3), S(i,4), S(i,5), S(i,6), c.l2]));
    v_g_squared = v_g_subbed(1)*v_g_subbed(1)+v_g_subbed(2)*v_g_subbed(2);
    integration_Kinetic_Energy = 1/2*c.m2*v_g_squared + 1/24*c.m2*c.l2^2*S(i,4)^2;
    integration_Potential_Energy = 1/2*c.k*(S(i,5)-c.L0)^2 - c.m2*c.g*(S(i,5)*cos(S(i,1))+c.l2/2*cos(S(i,3)));
    total_Energy = integration_Kinetic_Energy + integration_Potential_Energy;
    energy_vector4(i) = total_Energy;
end
energy_vector4 = energy_vector4(:,1);


figure(2)
hold on

grid on

(max(energy_vector3)-max(energy_vector4))/max(energy_vector4)*100

[min(energy_vector3), max(energy_vector3)]
ylim([min(energy_vector3), max(energy_vector3)])
plot(T,energy_vector3,'DisplayName', 'Default Tol')
plot(T,energy_vector4,'DisplayName', 'Increased Tol')
legend('show','Location','northoutside')
xlabel('Time, sec')
ylabel('Total Energy, J')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 6 2.5]);
print('-f2','Energy2','-dpdf');

hold off

%-------------------------------------------------------------------------%

T = T5;
S = S5;
energy_vector5 = zeros(length(S));
for i = 1:length(S)
    if mod(i,100)==0
        i
    end
    v_g_subbed = eval(subs(v_g, [thetat, thetadot, phit, phidot, l1t, l1dot, l2],...
        [S(i,1), S(i,2), S(i,3), S(i,4), S(i,5), S(i,6), c.l2]));
    v_g_squared = v_g_subbed(1)*v_g_subbed(1)+v_g_subbed(2)*v_g_subbed(2);
    integration_Kinetic_Energy = 1/2*c.m2*v_g_squared + 1/24*c.m2*c.l2^2*S(i,4)^2;
    integration_Potential_Energy = 1/2*c.k*(S(i,5)-c.L0)^2 - c.m2*c.g*(S(i,5)*cos(S(i,1))+c.l2/2*cos(S(i,3)));
    total_Energy = integration_Kinetic_Energy + integration_Potential_Energy;
    energy_vector5(i) = total_Energy;
end
energy_vector5 = energy_vector5(:,1);

T = T6;
S = S6;
energy_vector6 = zeros(length(S));
for i = 1:length(S)
    if mod(i,100)==0
        i
    end
    v_g_subbed = eval(subs(v_g, [thetat, thetadot, phit, phidot, l1t, l1dot, l2],...
        [S(i,1), S(i,2), S(i,3), S(i,4), S(i,5), S(i,6), c.l2]));
    v_g_squared = v_g_subbed(1)*v_g_subbed(1)+v_g_subbed(2)*v_g_subbed(2);
    integration_Kinetic_Energy = 1/2*c.m2*v_g_squared + 1/24*c.m2*c.l2^2*S(i,4)^2;
    integration_Potential_Energy = 1/2*c.k*(S(i,5)-c.L0)^2 - c.m2*c.g*(S(i,5)*cos(S(i,1))+c.l2/2*cos(S(i,3)));
    total_Energy = integration_Kinetic_Energy + integration_Potential_Energy;
    energy_vector6(i) = total_Energy;
end
energy_vector6 = energy_vector6(:,1);


figure(3)
hold on

grid on

(max(energy_vector5)-max(energy_vector6))/max(energy_vector6)*100

[min(energy_vector5), max(energy_vector5)]
ylim([min(energy_vector5), max(energy_vector5)])
plot(T,energy_vector5,'DisplayName', 'Default Tol')
plot(T,energy_vector6,'DisplayName', 'Increased Tol')
legend('show','Location','northoutside')
xlabel('Time, sec')
ylabel('Total Energy, J')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 6 2.5]);
print('-f3','Energy3','-dpdf');

hold off

%-------------------------------------------------------------------------%
