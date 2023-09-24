
function xbest=GROOVE_GA(B)
% for cont=1:1
% tic
info=[1, 0.025, 0.05, 21, 1; 
    2, 0, 1.39, 21, 1; 
    3, 0, 1.39, 21, 1; 
    4, 300, 1500, 21, 8; 
    5, 0.00025, 0.0005, 21, 1; 
    6, 0.17545, 1.0472, 21, 1; 
    7, 0.00025, 0.0005, 21, 1; 
    8, 0, 3, 2, 1; 
    9, 0, 1, 1, 1];
Np=100;         %Numero de particulas
K=50;           %Numero maximo de iteraciones
pc=0.3517;
Nb=21;
pm=0.5816;
k=0;
d=9;
%% INICIALIZACION: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B= 2.*rand(Np,(Nb*7)+3)-1;  %Binary population
% B=hardlim(B);
D=size(B,2);    %Size of candidate solution, in bits
X=zeros(Np,d);  %Real and integer population
f=zeros(1,Np);
%% EVALUACION DE POBLACION INICIAL: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ind=1:Np
    X(ind,:)=[DECOD(B(ind,1:21), info(1,2), info(1,3), info(1,4)), ...
        DECOD(B(ind,22:42),info(2,2),info(2,3),info(2,4)), ...
        DECOD(B(ind,(43:63)),info(3,2),info(3,3),info(3,4)), ...
        DECOD(B(ind,(64:84)),info(4,2),info(4,3),info(4,4)), ...
        DECOD(B(ind,(85:105)),info(5,2),info(5,3),info(5,4)), ...
        DECOD(B(ind,(106:126)),info(6,2),info(6,3),info(6,4)), ...
        DECOD(B(ind,(127:147)),info(7,2),info(7,3),info(7,4)), ...
        DECOD(B(ind,(148:149)),info(8,2),info(8,3),info(8,4)), ...
        DECOD(B(ind,150),info(9,2),info(9,3),info(9,4))];
    Co=funcion2(X(ind,:));
    if Co<9
        [TT, YY]=funcion(X(ind,:));
        [~, ind2]=max(TT);
        f(1,ind)=YY(ind2);
    else
        f(1,ind)=0.001;
    end
end
%% OBTENER OPTIMO GLOBAL:                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[f,ind]=sort(f, 'descend');
B=B(ind,:); 
X=X(ind,:);
fbest=f(1,1);
xbest=X(1,:);
bbest=B(1,:);
angulo=zeros(1,K);
while k<K
    %% SELECCION DE PADRES PARA CRUZA (METODO DE LA RULETA): %%%%%%%%%%%%%%
    E=sum(f); 
    E=f./E;
    q=zeros(1,Np);
    for ind=1:Np
        q(1,ind)=sum(E(1,1:ind));
    end
    padre1=zeros(1,Np); padre2=zeros(1,Np);
    for ind=1:Np
        padre1(1,ind)=RULETA(q);
        padre2(1,ind)=RULETA(q);
    end     
    %% CRUZA DE UN PUNTO (GENERAR HIJOS): %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Hb=[];
    for ind=1:Np
        hijo1=[]; 
        hijo2=[];
        r1=rand();
        if r1<=pc
            ptCruce=floor(r1*((Nb*7)+3)); 
            hijo1=[B(padre1(ind),1:ptCruce),B(padre2(ind),ptCruce+1:end)];
            hijo2=[B(padre2(ind),1:ptCruce),B(padre1(ind),ptCruce+1:end)];
        end
        Hb=[Hb;hijo1;hijo2];
    end     
    %% MUTACION: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sizeH=size(Hb,1);
    Hm=rand(sizeH,((Nb*7)+3))<pm;
    for ind=1:sizeH
        ind2=find(Hm(ind,:));
        Hb(ind,ind2)=~Hb(ind,ind2);
    end
    Hm=Hb;
    Hr=zeros(sizeH,d);
    for ind=1:sizeH
     Hr(ind,:)=[DECOD(Hm(ind,1:Nb), info(1,2), info(1,3), info(1,4)),...
     DECOD(Hm(ind,(Nb+1):(Nb*2)),info(2,2),info(2,3),info(2,4)),...
     DECOD(Hm(ind,(43:63)),info(3,2),info(3,3),info(3,4)),...
     DECOD(Hm(ind,((3*Nb)+1):(Nb*4)),info(4,2),info(4,3),info(4,4)),...
     DECOD(Hm(ind,((4*Nb)+1):(Nb*5)),info(5,2),info(5,3),info(5,4)),...
     DECOD(Hm(ind,((5*Nb)+1):(Nb*6)),info(6,2),info(6,3),info(6,4)),...
     DECOD(Hm(ind,((6*Nb)+1):(Nb*7)),info(7,2),info(7,3),info(7,4)),...
     DECOD(Hm(ind,((7*Nb)+1):((7*Nb)+2)),info(8,2),info(8,3),info(8,4)),...
     DECOD(Hm(ind,150),info(9,2),info(9,3),info(9,4))];
    end
    Hn=Hr;
    Hb=Hm;
    fHn=zeros(1,sizeH);
    for ind=1:sizeH
        Co=funcion2(Hn(ind,:));
        if Co<9
            [TT YY]=funcion(Hn(ind,:));
            [~,ind2]=max(TT);
            fHn(1,ind)=YY(ind2,1);
        else
            fHn(1,ind)=0.001;
        end
    end
    %% RESELECCION: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fHn=[f,fHn]; Hn=[X;Hn]; Hb=[B;Hb];
    [fHn,ind]=sort(fHn, 'descend');
    Hn=Hn(ind,:);
    Hb=Hb(ind,:);    
    if fHn(1,1)>fbest
        fbest=fHn(1,1);
        xbest=Hn(1,:);
        bbest=Hb(1,:);
    end
    k=k+1;
    angulo(k)=fbest;
    B=Hb(1:Np,:);
    f=fHn(1,1:Np);       
    X=Hn(1:Np,:);   
    disp(sprintf('generacion = %d, mejor fenotipo = %.6f', k, fbest))
end
xbest=[xbest,fbest];
%   A=toc;
% valores(cont,:)=opt1;
% mejores(cont)=fbest;
% mej(cont,:)=maxf;
% cont=cont+1
% end
% funcion1(xbest);
%figure
%plot(angulo)
end
