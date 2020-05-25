function rxn = ODHReactions(C_solid,Ts,R,Pt,Flowin,RxnKinetic,deltaS0,deltaH0,compnumber)
% This code is for modeling of ODH reaction kinetics

n_solid  = Flowin * C_solid ; % mole flow of each component [Nm^3/s * mol/m^3] = [mol/s]
nt_solid = sum(n_solid);

Ct_solid = sum(C_solid);   % total mole concentration in solid phase

% component order list: [C2H6 C2H4 O2 CO2 CO H2O N2]
P_solid = Pt*(C_solid/Ct_solid);

% component order list for reaction: [C2H6 C2H4 O2 CO2 CO H2O]
K = ones(1,6);
for n = 1:6
    K(n) = exp((deltaS0 - deltaH0*(1/Ts - 1/(25+273.15)) )/R );
end

Tetha_star = 1/(1 + K(1)*P_solid(1) + K(2)*P_solid(2) + sqrt(K(3)*P_solid(3)) + ...
                K(4)*P_solid(4) + K(5)*P_solid(5) + K(6)*P_solid(6) );
            
Tetha_O    = sqrt(K(3)*P_solid(3))/Tetha_star;
Tetha_C2H6 = K(1)*P_solid(1)/Tetha_star;
Tetha_C2H4 = K(2)*P_solid(2)/Tetha_star;

k = ones(1,5);
rxn = k(:);
for i = 1:5
    k(i)   = exp(RxnKinetic.Aprime(i) - RxnKinetic.EnergyA(i)/R *(1/Ts - 1/(25+273.15) ) );
    if i < 4    
        rxn(i) = k(i) * Tetha_O^(RxnKinetic.m) * Tetha_C2H6;
    else
        rxn(i) = k(i) * Tetha_O^(RxnKinetic.m) * Tetha_C2H4;
    end
end

% sum of rate of reactions
rxn = rxn * RxnKinetic.vcoffrxn(:,compnumber) ;

% sum of heat reactions

deltaH_reactants = nt_solid * Cpf * (298.15 - T0);
deltaH_products  = nt_solid * Cpf * (Ts - 298.15);
deltaH_std       = sum([RxnKinetic.deltaHstd]);

rxn = rxn * (20); 

end
