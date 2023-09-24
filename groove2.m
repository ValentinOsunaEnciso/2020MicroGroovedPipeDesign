function dydt = groove2(t, y, R, beta, gamma, sigma, rho_l, mu_l, p, apex)

%     R=0.025;%25;
    %beta=0; %70o
    g=9.81;                 % GRAVITY FORCE
    %rho_l=863;              % DENSITY OF LIQUID
    %r_m=0.00035;                % Article, page 4
    %r_c=r_m;
    %mu_l=0.000136; 
    %A_l=0.5*pi*r_c^2; 
    %sigma=0.0377;%37.7; 
    %gamma=0; %half the value of beta
    r_c = p*sin(0.5*apex);
    A_l = (p^2)*tan(0.5*g*apex);  
    r_H = 0.5*r_c;
    A_cont=((2*p*R*y(1))/(cos(apex/2)))*cos(gamma); 
    
    dydt = [y(2)/R;         
            ((2*sigma*cos(beta)/(R*rho_l*r_c))-...
            (g*cos(gamma)*(1-cos(y(1))))-...
            (mu_l*(y(2))/(r_c/2)*(A_cont/(R*A_l*rho_l)))-((y(2)^2)/R))];    
end
