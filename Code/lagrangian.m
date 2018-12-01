clear;close all;clc

% Constants
c.m2 = 1; % kg
c.k = 25; % N/m
c.L0 = 0.5;
c.l2 = 1; % m
c.g = 9.81; % m/s^2

syms l1 l2 m2 L0 g thetat thetadot thetaddot phit phidot phiddot k L
syms theta(t) phi(t)
a = [0 0 0];
g = [l1*sin(thetat)+l2/2*sin(phit) -l1*cos(thetat)-l2/2*cos(phit) 0];
r_ag = g - a;
v_g = subs(...
    diff(r_ag,thetat)*thetadot+diff(r_ag,phit)*phidot, [thetat, phit], [theta(t), phi(t)] ...
    );
v_g_squared = v_g(1)*v_g(1) + v_g(2)*v_g(2);

T1 = 0;
T2 = 1/2*m2*v_g_squared + 1/24*m2*l2*phidot^2;
T = T1 + T2;
V = 1/2*k*(l1-L)^2 - m2*g*(l1*cos(theta(t))+l2/2*cos(phi(t)));
L = T - V;
dLdthetadot = diff(L,thetadot)*thetaddot;
dLdthetadot_subbed = subs(dLdthetadot, [theta(t), phi(t)], [thetat, phit]);
ddLdthetadotdt = subs(...
    diff(dLdthetadot_subbed,thetat)*thetadot+diff(dLdthetadot_subbed,phit)*phidot, [thetat, phit], [theta(t), phi(t)] ...
    )
L_subbed = subs(L, [theta(t), phi(t)], [thetat, phit]);
dLdtheta = diff(L,thetat)*thetadot
