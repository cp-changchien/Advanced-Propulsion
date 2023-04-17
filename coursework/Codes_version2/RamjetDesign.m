function [A_1, A_C1, A_2, A_b, A_C2, A_4, eta_th, eta_p] = RamjetDesign(P_1, T_1, M_1, M_x, M_2, T_b, Pb_P2, P4_P1, F)
                
    % Function to calculate the cross-sectional areas at several key stations along a ramjet engine, 
    % together with its thermodynamic and propulsive efficiencies
    
    % INPUTS:   - P_1:    free-stream static pressure (Pa)
    %           - T_1:    free-stream temperature (K)
    %           - M_1:    flight Mach number
    %           - M_x:    normal shock strength
    %           - M_2:    burner entry Mach number
    %           - T_b:    burner temperature (K)
    %           - Pb_P2:  burner pressure ratio
    %           - P4_P1:  exhaust pressure ratio
    %           - F:      required thrust (N)
    
    % OUTPUTS:  - A_1:    inlet area (m^2)
    %           - A_C1:   inlet throat area (m^2)
    %           - A_2:    burner entry area (m^2)
    %           - A_b:    burner exit area (m^2)
    %           - A_C2:   nozzle throat area (m^2)
    %           - A_4:    exhaust area (m^2)
    %           - eta_th: thermodynamic efficiency
    %           - eta_p:  propulsive efficiency
    
    
    % STATIONS: - 1:      inlet (free-stream)
    %           - C1:     inlet throat
    %           - x:      immediately before shock
    %           - y:      immediately after shock
    %           - 2:      beginning of burner
    %           - b:      end of burner
    %           - C2:     nozzle throat
    %           - 4:      engine exhaust
    
    
    gamma = 1.4;  % Assume constant specific heat ratio throughout
    
    
    % INLET (1):
    [T01_T1, P01_P1, A1_AC1] = Isentropic(M_1, gamma);             % use isentropic relations to obtain temperature, pressure, and area ratios
    
    
    % SHOCK (x -> y):
    [Py_Px, M_y] = NormalShock(M_x, gamma);                        % use non-isentropic relations to obtain parameters across shock
    [~, P0x_Px, Ax_AxStar] = Isentropic(M_x, gamma);               % use isentropic relations to obtain pressure and area ratios
    [~, P0y_Py, Ay_AyStar] = Isentropic(M_y, gamma);               % use isentropic relations to obtain pressure and area ratios
    
    
    % BURNER ENTRY (2):
    [T02_T2, P02_P2, A2_A2Star] = Isentropic(M_2, gamma);          % use isentropic relations to obtain temperature, pressure, and area ratios
    A2_A1 = A2_A2Star * Ay_AyStar^(-1) * Ax_AxStar * A1_AC1^(-1);  % area ratio at burner entry, relative to inlet area
    T_2 = T02_T2^(-1) * T01_T1 * T_1;                              % static temperature at burner entry (note that T02 = T01 across shock)
    
    
    % ACROSS BURNER (2 -> b):
    % Equating conservation of mass and momentum across the burner, the following quadratic is found for burner exit Mach number (M_b):
    % M_b^2 - M_b * ( (M_2 + 1/(M_2*gamma))*sqrt(T_2/T_b) ) + 1/gamma = 0
    
    % Check for the sign of the discriminant b^2 - 4ac in the quadratic:
    discriminant = (T_2/T_b) * (M_2 + 1/(M_2*gamma))^2 - 4/gamma;
    if discriminant < 0
        disp("Invalid input, non-real burner exit Mach number");
        exit;
    elseif discriminant == 0
        M_b = 0.5*(M_2 + 1/(M_2*gamma))*sqrt(T_2/T_b);
    else
        M_b_solution1 = 0.5*(M_2 + 1/(M_2*gamma))*sqrt(T_2/T_b) + 0.5*sqrt(discriminant);  % 1st solution to the quadratic
        M_b_solution2 = 0.5*(M_2 + 1/(M_2*gamma))*sqrt(T_2/T_b) - 0.5*sqrt(discriminant);  % 2nd solution to the quadratic
        
        % Supersonic solution is impossible after the shock => choose subsonic solution (M_b < 1):
        if (M_b_solution1 < 1) && (M_b_solution2 > 1)
            M_b = M_b_solution1;
        elseif (M_b_solution1 > 1) && (M_b_solution2 < 1)
            M_b = M_b_solution2;
        else
            % If both solutions are subsonic, choose the weaker one:
            M_b = min([M_b_solution1, M_b_solution2]);
        end
    end
    
    % Finalise burner calculations:
    A2_Ab = Pb_P2 * (M_b/M_2) * sqrt(T_2/T_b);          % expression comes from mass conservation across the burner
    Ab_A1 = A2_Ab^(-1) * A2_A1;                         % area ratio at burner exit, relative to inlet area
    [T0b_Tb, P0b_Pb, Ab_AC2] = Isentropic(M_b, gamma);  % use isentropic relations to obtain temperature, pressure, and area ratios
    
    
    % NOZZLE THROAT (C2):
    AC2_A1 = Ab_AC2^(-1) * Ab_A1;                       % area ratio at nozzle throat, relative to inlet area
    
    
    % ENGINE EXHAUST (4):
    % Calculate stagnation-to-static pressure ratio at exit (noting P04/P0b = P02/P0y = P0x/P1 = 1 because of isentropic flow)
    P04_P4 = P0b_Pb * Pb_P2 .* P02_P2.^(-1) * P0y_Py * Py_Px .* P0x_Px.^(-1) * P01_P1 .* P4_P1.^(-1);
    
    % Finalise exhaust calculations:
    M_4 = sqrt( (2/(gamma-1)) * ((P04_P4).^((gamma-1)/gamma) - 1) );  % exhaust exit Mach number with isentropic relations
    [T04_T4, ~, A4_A4Star] = Isentropic(M_4, gamma);                  % use isentropic relations to obtain temperature and area ratios
    A4_A1 = A4_A4Star .* AC2_A1;                                      % area ratio at exhaust exit, relative to inlet area
    T_4 = T_b .* T0b_Tb .* T04_T4.^(-1);                              % static temperature at engine exhaust (note T_04/T_0b = 1)
    
    
    % THRUST:
    % From the momentum equation across the whole engine:
    F_A1 = gamma*P_1 * (P4_P1*A4_A1.*M_4.^2 - M_1.^2);  % F/A1, thrust normalised by inlet area
    
    
    % OUTPUT AREAS:
    A_1 = F ./ F_A1;           % inlet area
    A_C1 = A_1 * A1_AC1^(-1);  % inlet throat area
    A_2 = A2_A1 .* A_1;        % burner entry area
    A_b = Ab_A1 .* A_1;        % burner exit area
    A_C2 = AC2_A1 .* A_1;      % nozzle throat area
    A_4 = A4_A1 .* A_1;        % exhaust area
    
    
    % EFFICIENCIES:
    eta_th = 1 - (T_4 - T_1) / (T_b - T_2);                            % thermodynamic efficiency
    eta_p = ( 2*F_A1/(gamma*P_1) ) / ( (M_4.^2).*T_4/T_1 - (M_1.^2));  % propulsive efficiency
    
    
end