function [M2,rho2_rho1,U1_U2,P02_P01,P2_P1,T2_T1]=NormalShock(M1,gamma)
% function for calculating relations across a normal shock
% ------INPUTS------
% M1 - freestream mach number
% gamma - ratio of specific heats for air

%------------ Equations --------------
M2 = sqrt(((gamma-1)*M1^2+2)/(2*gamma*M1^2-(gamma-1))); % Mach number downstream of the shock.
rho2_rho1 = ((gamma+1)*M1^2)/((gamma-1)*M1^2+2); % density change across the shock.
U1_U2 = rho2_rho1; % velocity change across the shock
P02_P01 = ((2*gamma*M1^2-(gamma-1))/(gamma+1))^(-1/(gamma-1))*(((gamma+1)*M1^2)/(2+(gamma-1)*M1^2))^(gamma/(gamma-1)); % stagnation pressure change across the shock
P2_P1 = (2*gamma*M1^2-(gamma-1))/(gamma+1); % static pressure change across the shock
T2_T1 = ((2*gamma*M1^2-(gamma-1))*(2+(gamma-1)*M1^2))/((gamma+1)^2*M1^2); % temperature change across the shock