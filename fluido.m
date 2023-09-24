% Valentín Osuna-Enciso, Mayo, Junio, 2016; CUTONALA-UDG.
% Eqs. DIPPR105, DIPPR106, Vogel Equation
function [rho,sigma,mu] = fluido(kf,T)
%% T: temperature, oC, 150o <=T<=239o
%  kf: kind of fluid: 1)benzene, 2)water, 3)ethanol, 4)methanol
%  returns:
%  rho: density, in kg/m^3
%  sigma: surface tension, in N/m
%  mu: dinamic viscosity of liquid, in Pa*s
    switch (kf)
        case 1 % Benzene         
            A1=80.3682; B1=0.266457; C1=561.65; D1=0.28797;            
            A2=72.989; 	 B2=1.2578; 	 C2=0; 	 D2=0; 	 E2=0; 	 Tc=562.1;            
            A3=-7.92779;	B3=4313.54;	C3=283.29;            
        case 2 % Water
            A1=0.14395; 	 B1=0.0112; 	 C1=649.727; 	 D1=0.05107;            
            A2=134.15;  	 B2=1.6146;  	 C2=-2.035;  	 D2=1.5598;  	 E2=0;  	 Tc=647.3;            
            A3=-3.7188;	B3=578.919;	C3=-137.546;	
        case 3 % Ethanol 
            A1=99.3974; 	 B1=0.310729; 	 C1=513.18; 	 D1=0.305143;            
            A2=131.38; 	 B2=5.5437; 	 C2=-8.4826; 	 D2=4.3164; 	 E2=0; 	 Tc=516.2;            
            A3=-7.37146;	B3=2770.25;	C3=74.6787;	            
        case 4 %Methanol
            A1=54.566; 	 B1=0.233211; 	 C1=513.16; 	 D1=0.208875;            
            A2=102.06; 	 B2=4.2709; 	 C2=-6.0509; 	 D2=2.9715; 	 E2=0; 	 Tc=512.6;            
            A3=-6.7562;	B3=2337.24;	C3=84.0853;	            
        otherwise
            disp('Unknown liquid.')  
    end
    rho=A1/(B1^(1+(1-((T+273.15)/C1))^D1));
    sigma=(A2/1000)*(1-((T+273.15)/Tc))^(B2+(C2*((T+273.15)/Tc))+(D2*((T+273.15)/Tc)^2)+(E2*((T+273.15)/Tc)^3));
    mu=(exp(A3+(B3/(C3+(T+273.15)))))/1000;
end


