%% DECODIFICA CADENAS DE BITS %%
% Tamano de paso: ((u - l) / t3)
function r = DECOD1(B,l,u,Nb)
    % r:    valor real
    % B:    cadena binaria de Nb bits
    % Nb:   numero de bits de la cadena binaria
    % l:    limite inferior del espacio de busqueda
    % u:    limite superior del espacio de busqueda 
    r=zeros(1,size(l,2)+2); 
    i2=1; i3=Nb;
    for i1=1:size(l,2)
        temp = fliplr(B(i2:i3));        
        i2=i2+Nb; i3=i3+Nb;
        t1=(Nb-1):-1:0;
        t2=ones(1,Nb);          % Numero binario maximo a obtener
        t3=sum((t2.*2.^t1)')+1; % Numero de posibles valores reales
        t2=sum((temp.*2.^t1)');    % Numero real del individuo
        r(1,i1) =  l(1,i1) + t2 .* ((u(1,i1) - l(1,i1)) / t3);  % Aplicar limites del espacio de busque
    end
    temp = B(end-2:end-1); i1=i1+1;
    if temp(1,1)==0 && temp(1,2)==0
        r(1,i1)=1;
    elseif temp(1,1)==0 && temp(1,2)==1
        r(1,i1)=2;
    elseif temp(1,1)==1 && temp(1,2)==0
        r(1,i1)=3; 
    else
        r(1,i1)=4; 
    end
    temp = B(end); i1=i1+1;
    if temp==0
        r(1,i1)=1;
    else
        r(1,i1)=2;
    end
end