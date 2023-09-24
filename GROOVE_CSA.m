% Clonal Selection Algorithm: Design of microgrooved pipe in a parabolic 
% trough, CUTONALA, Jos� Valent�n Osuna-Enciso, Mayo, Junio, 2016.
function A=GROOVE_CSA(B)   
    tic
    %% PAR�METROS:    
    Np = 100;       % Tama�o de la poblaci�n
    n = 90;         % N�mero de individuos a clonar
    L=150;           % N�mero total de bits en el individuo
    Nb = 21;        % Bits por cada dimensi�n del individuo  
    maxk = 50;      % N�mero m�ximo de iteraciones
    pm = 0.5583;      % Probabilidad de mutaci�n
    Pd = 0.4383;      % Porcentaje de diversidad
    Pc = 0.1;      % Porcentaje de clonaci�n             
    k = 0;          % Iteraci�n actual
    nc = Np*Pc*n;   % N�mero total de clones en cada iteraci�n
    l = [0.025, 0,     0,     200, 0.00025, 0.1745, 0.00025];
    u = [0.05, 1.39, 1.39, 330, 0.00043, 1.0472, 0.00043];
    %% INICIALIZACI�N:
%     B = 2 .* rand(Np,L) - 1;   
%     B = hardlim(B);                  
    X=[]; fx = [];                       
    for i1 = 1 : Np
       X=[X;DECOD1(B(i1,:),l,u,Nb)];
       fx = [fx; rungroove(X(i1,:));];
    end
    ind=1:Np;
    %% ALGORITMO DE SELECCI�N CLONAL:
    while k <=  maxk
        %valx = X(ind(1:end),:); 
        %% CLONACI�N:
        [C,pcs] = CLONES(n,Pc,Np,ind,B);
        %% HIPERMUTACI�N:
        Cm=C;
        for i1=1:nc
            for i2=1:L
                if rand()<=pm
                    Cm(i1,i2)=~C(i1,i2);              
                end
            end
        end
        Cm(pcs,:) = B(Np-n+1:end,:);%B(ind(Np-n+1:end),:);              
        Xm=[]; fm = [];                                 
        for i1 = 1 : nc
            Xm=[Xm;DECOD1(Cm(i1,:),l,u,Nb)];
            fm = [fm; rungroove(Xm(i1,:));];            
        end    
        %% RESELECCI�N:
        pcs = [0 pcs];
        for i1=1:n
            [out(i1),bcs(i1)] = min(fm(pcs(i1)+1:pcs(i1+1)));  
            bcs(i1) = bcs(i1) + pcs(i1);
        end
        B(Np-n+1:end,:) = Cm(bcs,:);
        X(Np-n+1:end,:) = Xm(bcs,:);
        fx(Np-n+1:end,1)=out;
        [fx,ind] = sort(fx);
        temp=X(ind(1:end),:); X=temp;
        temp=B(ind(1:end),:); B=temp;
        %% INTRODUCCI�N DE DIVERSIDAD:   
        nd = round(Pd*Np); 
        B(ind(end-nd+1:end),:) = 2 .* rand(nd,L) - 1;   
        B(ind(end-nd+1:end),:) = hardlim(B(ind(end-nd+1:end),:)); 
        %% MUESTRA AVANCES EN PANTALLA:
        k = k + 1;
        %pm = pmcont(pm,pma,pmr,k,maxk); 
        valfx(k) =fx(1);
        fbest=fx(1);
        valx=X(1,:);
        disp(sprintf('It:%d,R:%2.2f,b:%2.2f,g:%2.2f,st:%2.2f,p:%2.2f,psi:%2.2f,rm:%2.2f,kf:%d,tg:%d,f:%2.5f',k,valx,fbest));        
        %rungroove_imp(valx)
    end % end while
    fbest=rungrooved(valx);
    A=toc;
%     plot(valfx), 
%     rungroove_imp(valx)
end
