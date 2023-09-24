% Binary Differential Evolution, 
function xbest=GROOVE_DBDE(B)
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
    F=0.5025;  %Step size,         [0, 2]
    CR=0.5431; %Crossover rate,    [0, 1]
    Np=100; %Population size
    K=50;   %Maximum iteration number
    k=0;    %Actual iteration
    Nb=21;  %Number of bits per variable in candidate solution  
    angulo=zeros(1,K);
%     B= 2.*rand(Np,(Nb*7)+3)-1;   %Binary population
%     B=hardlim(B);
    D=size(B,2);    %Size of candidate solution, in bits
    X=zeros(Np,9);  %Real and integer population
    f=zeros(1,Np);
    for ind=1:Np
        X(ind,:)=[DECOD(B(ind,1:21), info(1,2), info(1,3), info(1,4)),...
            DECOD(B(ind,22:42),info(2,2),info(2,3),info(2,4)), ...
            DECOD(B(ind,(43:63)),info(3,2),info(3,3),info(3,4)),...
            DECOD(B(ind,(64:84)),info(4,2),info(4,3),info(4,4)),...
            DECOD(B(ind,(85:105)),info(5,2),info(5,3),info(5,4)),...
            DECOD(B(ind,(106:126)),info(6,2),info(6,3),info(6,4)),...
            DECOD(B(ind,(127:147)),info(7,2),info(7,3),info(7,4)),...
            DECOD(B(ind,(148:149)),info(8,2),info(8,3),info(8,4)),...
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
    [fbest,ind]=max(f); xbest=X(ind,:);
    v=zeros(1,D); u=zeros(1,D);
    while K>k
        for jj=1:Np            
            r=randperm(Np);
            jrand=randi([1,D]); 
            r1=r(1);
            r2=r(2);
            r3=r(3);
            xr1=B(r1,1:150);
            xr2=B(r2,1:150);
            xr3=B(r3,1:150);
            % Mutation:
            v=xr3+F.*(xr1-xr2);
            for d=1:D
                 if rand < (1/(1+exp(-v(1,d))))
                     v(1,d)=1;
                 else
                     v(1,d)=0;
                 end
            end                            
            % Crossover:
            for d=1:D
                if d==jrand || rand() <= CR
                    u(1,d)=v(1,d);
                else
                    u(1,d)=B(jj,d);
                end 
            end
            U=[DECOD(u(1,1:21), info(1,2), info(1,3), info(1,4)), ...
                DECOD(u(1,22:42),info(2,2),info(2,3),info(2,4)), ...
                DECOD(u(1,(43:63)),info(3,2),info(3,3),info(3,4)),...
                DECOD(u(1,(64:84)),info(4,2),info(4,3),info(4,4)),...
                DECOD(u(1,(85:105)),info(5,2),info(5,3),info(5,4)),...
                DECOD(u(1,(106:126)),info(6,2),info(6,3),info(6,4)),...
                DECOD(u(1,(127:147)),info(7,2),info(7,3),info(7,4)),...
                DECOD(u(1,(148:149)),info(8,2),info(8,3),info(8,4)),...
                DECOD(u(1,150),info(9,2),info(9,3),info(9,4))];             
            % Selection:
            Co=funcion2(U);
            if Co<9
                [TT,YY]=funcion(U);
                [~, ind2]=max(TT);
                fU=YY(ind2);
            else
                fU=0.001;
            end
            if fU > f(jj)
               f(1,jj)=fU;
               B(jj,:)=u;
               X(jj,:)=U;
            end          
            if fU>fbest
               xbest=U;
               fbest=fU;
            end
        end
        k=k+1;  
        angulo(k)=fbest;
        disp(sprintf('generacion = %d, mejor fenotipo = %.6f', k, fbest))
    end
    xbest=[xbest,fbest];
%     A=toc;
%     disp(['Mejor posición x:',num2str(xbest)])
%     funcion1(xbest);
%     figure
%     plot(angulo);
end