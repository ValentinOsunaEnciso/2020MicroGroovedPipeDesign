function r = DECOD(B,l,u,Nb)
%r: valor real
%B Cadena binaria de Nb bits
%Nb n�mero de bits de la cadena binaria
%l l�mite inferior del espacio de b�squeda
%u l�mite superior del espacio de b�squeda
if Nb<2
    r=bi2de(B);
elseif    Nb<3
    r=bi2de(B);
else
B=fliplr(B);
t1=(Nb-1):-1:0;
t2=ones(1,Nb); %N�mero binario m�ximo a obtener
t3=sum((t2.*2.^t1)')+1; %n�mero de posibles valores reales
t2=sum((B.*2.^t1)');   %N�mero real del individuo
r = l+t2.*((u-l)/t3); %Aplicar l�mites del espacio de busqueda
end
end