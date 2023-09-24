function r = DECOD(B,l,u,Nb)
%r: valor real
%B Cadena binaria de Nb bits
%Nb número de bits de la cadena binaria
%l límite inferior del espacio de búsqueda
%u límite superior del espacio de búsqueda
if Nb<2
    r=bi2de(B);
elseif    Nb<3
    r=bi2de(B);
else
B=fliplr(B);
t1=(Nb-1):-1:0;
t2=ones(1,Nb); %Número binario máximo a obtener
t3=sum((t2.*2.^t1)')+1; %número de posibles valores reales
t2=sum((B.*2.^t1)');   %Número real del individuo
r = l+t2.*((u-l)/t3); %Aplicar límites del espacio de busqueda
end
end