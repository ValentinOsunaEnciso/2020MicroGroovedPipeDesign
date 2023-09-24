
function [TT U]= funcion(s)

 
R=s(1);
B=s(2);
v=s(3); 
st=s(4);
p=s(5);
ang=s(6);
rm=s(7);
kf=s(8);
tg=s(9);
 %st  Temperatura de saturacion grados centigrados
 %R   Radio interno del tubo
 %o   tensión de la superficie
 %B   Angulo de contacto de la interfaz
 %d;  Densidad del liquido
 %p;  altura del surco
 %rc; radio efectivo capilar equivalente en la dirección de interes
 %m  Viscocidad dinamica del liquido
 g = 9.81; %aceleración de la gravedad
if kf==0                    %LIQUID SODIUM
    if st<371
       st=371; 
    end
    if st>1600
       st=1600; 
    end
    d=219+(275.32*(1-(st/2503.7)))+(511.58*((1-(st/2503.7))^0.5));
    o=240.5*((1-(st/2503.7))^1.126);
    m=(exp(-6.4406-(0.3958*log(st))+(556.835/st)))/1000;
elseif kf==1                %WATER
    A1= 0.14395;
    B1= 0.0112;
    C1= 649.727;
    D1= 0.05107;
    A2= 134.15;
    B2= 1.6146;
    C2= -2.035;
    D2= 1.5598;
    E2= 0;
    Tc= 647.3;
    A3= -3.7188;
    B3= 578.919;
    C3= -137.546;
    if st<233
       st=233; 
    end
    if st>643
       st=643; 
    end
    d=A1/(B1^(1+((1-(st/C1))^D1)));
    o=(A2/1000)*((1-(st/Tc))^(B2+(C2*(st/Tc))+(D2*((st/Tc)^2))+(E2*((st/Tc)^3))));
    m=(exp(A3+(B3/(C3+st))))/1000;
  elseif kf==2              %HITEC SALT
    if st<450
       st=450; 
    end
    if st>1050
       st=1050; 
    end  
    d=2293.6-(0.7497*st);
    o=(0.14928-((0.556^-4)*st))/1000;
    m=(0.4737-(0.002297*st)+(3.3731^-6*st^2)-(2.019^-9*st^3))/1000;
elseif kf==3                %MOLTEN SALT
    if st<750
       st=750; 
    end
    if st>1550
       st=1550; 
    end
    d=2363.84-(0.474*st);
    o=(0.133-(0.48^-4*st))/1000;
    m=((1.46^-4)*exp(2230/st))/1000;
end  

[TT,YY]=call_osc();
function dydt = osc(t,y)
    
if tg==0   
    Acont = (pi*rm*R*y(1)*cos(B));
    A_1 = (1/2)*pi*rm^2;
    rc=rm;
    rh= (1/2)*rm;
 else
    Acont = ((2*p*R*y(1))/(cos(ang/2)))*cos(B);
    A_1 = (p^2)*(tan(ang/2));
    rc = p*sin(ang/2);
    rh= (1/2)*p*sin(ang/2);
 end

 dydt = zeros(2,1); % this creates an empty column
 %vector that you can fill with your two derivatives:
 dydt(1) = y(2)/R;
 dydt(2) = ((2*o*cos(B))/(R*d*rc))-(g*cos(v)*(1-cos(y(1))))-(m*(y(2)/(rh))*(Acont/(R*A_1*d)))-((y(2))^2/R);
 %In this case, y(1) is y1 and y(2) is y2, and dydt(1)
 %is dy1/dt and dydt(2) is dy2/dt.

end

function [T,Y] = call_osc()
 tspan = [0 10];
 y1_0 = 0;
 y2_0 = 0;
 [T,Y] = ode15s(@osc,tspan,[y1_0 y2_0]); %15s original
end
U=YY(:,1);
end

