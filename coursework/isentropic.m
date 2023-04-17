function [T01_T1, P01_P1, A_Achoked] = isentropic(M_1, gamma)
% isentropic relation for ranjet design 
% required input:
% M_1, Flight/operation Mach number
% gamma, gas constant (assumed to be 1.4)

T01_T1=1+(M_1^2)*(gamma-1)/2;
P01_P1=(1+(M_1^2)*(gamma-1)/2)^(gamma/(gamma-1));
A_Achoked=(1/M_1)*(((2/(gamma+1))*(1+(M_1^2)*(gamma-1)/2))^((gamma+1)/(2*(gamma-1))));

end