
corridas=200;
resultados=zeros(corridas,30); 
for ind=1:corridas                  %numero de corridas
Np=100;
Nb=21;
B= 2.*rand(Np,(Nb*7)+3)-1;      %Binary population
B=hardlim(B);
resultados(ind,:)=[GROOVE_GA(B),GROOVE_BPSO(B),GROOVE_DBDE(B)];
SALIDA=resultados(ind,:);
save('EXPERIMENTO_1.txt','SALIDA','-ASCII','-append');
end
Np=100;