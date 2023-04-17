function [T0_T, P0_P, A_AStar] = Isentropic(M, gamma)

    % Function to calculate isentropic flow relations
    
    % INPUTS:  - M:       Mach number at which to evaluate the relations
    %          - gamma:   specific heat ratio
    
    % OUTPUTS: - T0_T:    ratio of stagnation to static temperature at input Mach number 
    %          - P0_P:    ratio of stagnation to static pressure at input Mach number
    %          - A_AStar: ratio of arbitrary area (where Mach = M) to choked area (where Mach = 1)
    
    
    % Calculations:
    T0_T = 1 + 0.5*(gamma-1).*M.^2;
    P0_P = (1 + 0.5*(gamma-1).*M.^2) .^ (gamma/(gamma-1));
    A_AStar = (1./M).*( (2./(gamma+1)).*(1 + 0.5*(gamma-1).*M.^2) ) .^ (0.5*(gamma+1)/(gamma-1));

end