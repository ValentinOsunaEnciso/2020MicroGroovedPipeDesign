
corridas=30;
resultados=zeros(corridas,4);    
for ind=1:corridas                  %numero de corridas
Np=100;
Nb=21;
B= 2.*rand(Np,(Nb*7)+3)-1;      %Binary population
B=hardlim(B);
resultados(ind,:)=[GROOVE_CSA(B),GROOVE_GA(B),GROOVE_BPSO(B),GROOVE_DBDE(B)];
end
Np=100;