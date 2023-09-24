function dydt = groove1(t, y, R, beta, gamma, sigma, rho_l, r_m, mu_l)

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
    r_c = r_m;
    A_l = 0.5*pi*r_c^2;   
    r_H = 0.5*r_c;
    A_cont=pi*r_c*R*y(1)*cos(sigma); 
    
    dydt = [y(2)/R;         
            ((2*sigma*cos(beta)/(R*rho_l*r_c))-...
            (g*cos(gamma)*(1-cos(y(1))))-...
            (mu_l*(y(2))/(r_c/2)*(A_cont/(R*A_l*rho_l)))-((y(2)^2)/R))];    
end
