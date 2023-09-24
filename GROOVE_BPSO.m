function xbest=GROOVE_BPSO(B)    
% tic
%% PARAMETROS INICIALES: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info=[1, 0.025, 0.05, 21, 1; 
    2, 0, 1.39, 21, 1; 
    3, 0, 1.39, 21, 1; 
    4, 300, 1500, 21, 8; 
    5, 0.00025, 0.0005, 21, 1; 
    6, 0.17545, 1.0472, 21, 1; 
    7, 0.00025, 0.0005, 21, 1; 
    8, 0, 3, 2, 1; 
    9, 0, 1, 1, 1];
d=9;            %Size of candidate solution, in real numbers
c1=0.9020;
c2=0.5425;
w=1.2175;
Np=100;         %Population size
Nb=21;          %Number of bits per variable in candidate solution 
K=50;           %Maximum iteration number
v1=zeros(Np,150);
v0=zeros(Np,150);
angulo=zeros(1,K);
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
%% OBTENER OPTIMO GLOBAL Y OPTIMO LOCAL: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fbest,ind]=max(f);
bbest=B(ind,:);
xbest=X(ind,:);
blocal=B;
flocal=f;
%% ITERACIONES: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v=zeros(Np,150);
k=0;        %Actual iteration
while k<K
  for i=1:Np
    btemp=zeros(1,D);  
    for j=1:D
        r1=rand();
        r2=rand();
        if blocal(i,j)==1
            d11=c1*r1;
            d01=-c1*r1;
        else
            d01=c1*r1;
            d11=-c1*r1;
        end
        if bbest(1,j)==1
            d12=c2*r2;
            d02=-c2*r2;
        else
            d02=c2*r2;
            d12=-c2*r2;
        end                
        v1(i,j)=w*v1(i,j)+d11+d12;  %Equation (6)
        v0(i,j)=w*v0(i,j)+d01+d02;  %Equation (7)
        if B(i,j)==0
            vc=v1(i,j);
            r=rand();
            if r < (1/(1+exp(-vc)))
                btemp(1,j)=1;
            end
        else 
            vc=v0(i,j);
            r=rand();
            if r < (1/(1+exp(-vc)))
                btemp(1,j)=0;
            end
        end
    end    
    tempX=[DECOD(btemp(1,1:Nb), info(1,2), info(1,3), info(1,4)),...
       DECOD(btemp(1,(Nb+1):(Nb*2)),info(2,2),info(2,3),info(2,4)),...
       DECOD(btemp(1,(43:63)),info(3,2),info(3,3),info(3,4)),...
       DECOD(btemp(1,((3*Nb)+1):(Nb*4)),info(4,2),info(4,3),info(4,4)),...
       DECOD(btemp(1,((4*Nb)+1):(Nb*5)),info(5,2),info(5,3),info(5,4)),...
       DECOD(btemp(1,((5*Nb)+1):(Nb*6)),info(6,2),info(6,3),info(6,4)),...
       DECOD(btemp(1,((6*Nb)+1):(Nb*7)),info(7,2),info(7,3),info(7,4)),...
       DECOD(btemp(1,((7*Nb)+1):((7*Nb)+2)),info(8,2),info(8,3),info(8,4)),...
       DECOD(btemp(1,150),info(9,2),info(9,3),info(9,4))];
    Co=funcion2(tempX);   
    if Co<9
        [TT, YY]=funcion(tempX);       
        [~, ind2]=max(TT);
        tempf=YY(ind2);
    else
        tempf=0.001;
    end  
    %% COMPARAR PERFORMANCES: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if tempf>flocal(1,i)
        flocal(1,i)=tempf;
        blocal(i,:)=btemp;
        B(i,:)=btemp;
        X(i,:)=tempX;
        f=tempf;
        if tempf>fbest
            fbest=tempf;
            bbest=btemp;
            xbest=tempX;
        end
    end               
  end
  k=k+1;
  angulo(k)=fbest;
  disp(sprintf('generacion = %d, mejor fenotipo = %.6f', k, fbest))
end 
xbest=[xbest,fbest];
%   A=toc;
% funcion1(xbest);
% figure
% plot(angulo);
% title('Convergencia')
% ylabel('angulo')
% xlabel('iteracion')
end
