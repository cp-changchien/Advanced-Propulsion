function [Py_Px, M_y] = NormalShock(M_x, gamma)

    % Function to calculate non-isentropic flow relations across a normal shock
    
    % INPUTS:  - M_x:   Mach number immediately before the shock
    %          - gamma: specific heat ratio
    
    % OUTPUTS: - Py_Px: ratio of static pressures across the shock
    %          - M_y:   Mach number immediately after the shock
    
    
    % Calculations:
    Py_Px = (2*gamma*M_x.^2 - (gamma-1)) / (gamma+1);
    M_y = sqrt( ((gamma-1)*M_x.^2 + 2) / (2*gamma*M_x.^2 - (gamma-1)) );

end