function [A1, At1, Abe, Abx, At3, Ae, thermeff, propeff] = ramjet(pa,Ta,Mf,Shock,Mbe,Tb,bpr,epr,thrust)
        
        %----- Pre-defined values to test --------
        % pa = 7e4; %ambient pressure
        % Ta = 210; %ambient temperature
        % Mf = 2.8; %free stream mach number
        % Shock = 1.1; %normal shock strength (before shock mach number)
        % Mbe = 0.2; %burner entry mach number
        % Tb = 1700; %maximum burner temperature
        % bpr = 1; %burner pressure ratio
        % epr = 1; %exhaust pressure ratio
        % thrust = 5e5; %required thrust

        %------ Explanation on index -------------
        %a - ambient, 1 - inlet, 01/t1 - stagnation/choke conds before shock, 2b - just
        %before shock, 2a - just after shock, 02/t2 - stagnation/choke conds after shock,
        %before burner, be - burner entry, bx - burner exit, 03/t3 -
        %stagnation/choke conds after burner, e - exhaust
        %-----------------------------------------
        
        gamma = 1.4;

        [T01Ta,p01pa,A1At1] = isentropic(Mf,gamma);                 %inlet stagnation conditons and inlet throat area
        [T01T2b,p01p2b,A2bAt1] = isentropic(Shock,gamma);           %stagnation conditions before shock same as inlet (isentropic)
        A2bA1 = A2bAt1/A1At1;                                       %area ratio at shock
        T01 = T01Ta * Ta;                                           %stag temp before shock
        T02 = T01;                                                  %shock doesnt change stag temp
        T2b = T01/T01T2b;                                           %temp before shock
        p01 = pa*p01pa;                                             %stag pressure before shock
        p2b = p01/p01p2b;                                           %pressure just before shock
        [M2a,~,~,p02p01,p2ap2b,T2aT2b] = NormalShock(Shock,gamma);  %calculate all shock ratios
        p02 = p02p01*p01;                                           %stag pressure after shock
        p2a = p2ap2b*p2b;                                           %pressure just after shock
        T2a = T2b*T2aT2b;                                           %temp just after shock
        [T02T2a, p02p2a, A2aAt2] = isentropic(M2a,gamma);           %isentropic until burner entry
        [T02Tbe, p02pbe, AbeAt2] = isentropic(Mbe,gamma);


        %burner entry conditions
        AbeA2a = AbeAt2/A2aAt2;
        AbeA1 = AbeA2a*A2bA1;
        Tbe = T2a*T02T2a/T02Tbe;
        pbe = p2a*p02p2a/p02pbe;


        %burner exit conditions
        pbx = pbe*bpr;
        Tbx = Tb;


        %calculate burner exit mach number
        %check for imaginary root
        root = Tbe/Tbx * (Mbe + 1/gamma/Mbe)^2 - 4/gamma;
        if root<=0
            disp("Invalid input, solution unavailable");
            exit;
        end


        %choose subsonic solution
        Mbx1 = 0.5*sqrt(Tbe/Tbx) * (Mbe + 1/gamma/Mbe) + 0.5*sqrt(root);
        Mbx2 = 0.5*sqrt(Tbe/Tbx) * (Mbe + 1/gamma/Mbe) - 0.5*sqrt(root);
        if Mbx1 < 1 && Mbx1 > 0
            Mbx = Mbx1;
        elseif Mbx2 < 1 && Mbx2 > 0
            Mbx = Mbx2;
        else
            disp("Invalid input, solution unavailable");
        end


        %burner exit area ratios
        AbxAbe = (gamma*Mbe^2 + 1)/(bpr*(gamma*Mbx^2 + 1));
        AbxA1 = AbeA1*AbxAbe;


        %isentropic until exhaust
        %calculate stagnation and throat conditions
        [T03Tbx, p03pbx, AbxAt3] = isentropic(Mbx,gamma);
        p03 = p03pbx*pbx;
        At3A1 = AbxA1/AbxAt3;
        pe = epr*pa;
        p03pe = p03/pe;


        %solve for exhaust mach number from isentropic relation
        syms mach
        Mesolve = double(solve(p03pe == (1+(mach^2)*(gamma-1)/2)^(gamma/(gamma-1)),mach));
        

        %only choose real positive value
        for i = 1:length(Mesolve)
            if real(Mesolve(i)) > 0 && imag(Mesolve(i)) == 0
                Me = Mesolve(i);
            end
        end
        
        %isentropic relations at exhaust - same stagnation/choke as burner exit
        [T03Te, ~, AeAt3] = isentropic(Me,gamma);
        Te = T03Tbx/T03Te * Tbx; %exhaust temp
        AeA1 = AeAt3*At3A1; %area ratio


        %solve for inlet area from thrust required
        syms A1sym
        A1 = double(solve(thrust/(pa*A1sym) == gamma*Mf^2*(Me^2*AeA1/Mf^2 - 1)));
        
        
        %calculate all other outputs
        At1 = A1/A1At1;
        Abe = AbeA1*A1;
        Abx = AbxA1*A1;
        At3 = At3A1*A1;
        Ae = AeA1*A1;

        %efficiencies
        thermeff = 1 - (Te - Ta)/(Tbx - Tbe);
        propeff = thrust/(pa*A1) * 2*Ta/(gamma*(Te*Me^2 - Ta*Mf^2));
        
end