v_opt=2.61;
EXPERIMENTOS=zeros(200,10);
for ind3=1:200
info=[1, 0.025, 0.05, 21, 1; 
    2, 0, 1.39, 21, 1; 
    3, 0, 1.39, 21, 1; 
    4, 421, 501, 21, 8; 
    5, 0.00025, 0.0005, 21, 1; 
    6, 0.17545, 1.0472, 21, 1; 
    7, 0.00025, 0.0005, 21, 1; 
    8, 0, 3, 2, 1; 
    9, 0, 1, 1, 1];
L=[0;1.08016760892075e-06;3.40913290388060e-06;7.01321568223160e-06;3.66637039233811e-05;9.01475160263803e-05;0.000167838439240867;0.000269865509733788;0.000721420418964088;0.00139145137819088;0.00227966792607881;0.00338567057045700;0.00816521283141612;0.0149914716046106;0.0238348026997321;0.0346564239960897;0.0586188325563506;0.0883558148798499;0.123467396813599;0.163490857945160;0.207919695912044;0.303584189525487;0.409728802505654;0.522698322220427;0.639203160770772;0.756480683578235;0.943598082569695;1.12053514772076;1.28339497835791;1.43054347170903;1.56177642756835;1.68616506745487;1.79410351252050;1.88729762092726;1.96755899690895;2.03664943306928;2.09617988886447;2.16542947137854;2.22224909535721;2.26928097351234;2.30858721514515;2.34175590873659;2.37449909880282;2.40197905345372;2.42529607580304;2.44528404898472;2.46257873879393;2.47766265809052;2.49533929546870;2.51028759745192;2.52302600627890;2.53395483237409;2.54338509929041;2.55155886904146;2.56094190015609;2.56875274089676;2.57528656384892;2.58077394630154;2.58539527500667;2.58929551630014;2.59259400569108;2.59602791013500;2.59882453828184;2.60110526009428;2.60296772624958;2.60449078718713;2.60573769857029;2.60737783990090;2.60854570111179;2.60937645836990;2.60996737710448;2.61038842203070;2.61068903186386;2.61094078263094;2.61110841801593;2.61121945243634;2.61129297200820;2.61134195058596;2.61137472765995;2.61141472496374;2.61143234602844;2.61143787298144;2.61143848324331;2.61143866201615;2.61143941234888;2.61144001573882;2.61144010305062;2.61143992189338;2.61143980609097;2.61143981672285;2.61143982604784;2.61143985861968;2.61143984949268;2.61143983764590;2.61143983620825;2.61143983762664];
Np=100;
Nb=21;   
k=1;
K=51;
%Inicialización
B= 2.*rand(Np,(Nb*7)+3)-1;  %Binary population
B=hardlim(B);
T=[];
R=[];
for ind=1:Np
    R=[R;DECOD(B(ind,1:21), info(1,2), info(1,3), info(1,4)), DECOD(B(ind,22:42),info(2,2),info(2,3),info(2,4)), DECOD(B(ind,(43:63)),info(3,2),info(3,3),info(3,4)),DECOD(B(ind,(64:84)),info(4,2),info(4,3),info(4,4)),DECOD(B(ind,(85:105)),info(5,2),info(5,3),info(5,4)),DECOD(B(ind,(106:126)),info(6,2),info(6,3),info(6,4)),DECOD(B(ind,(127:147)),info(7,2),info(7,3),info(7,4)),DECOD(B(ind,(148:149)),info(8,2),info(8,3),info(8,4)),DECOD(B(ind,150),info(9,2),info(9,3),info(9,4))];
end 
%%%%% Evaluaciòn de la funciòn
for ind1=1:Np
     Co=funcion2(R(ind1,:));
     if Co<9
     [TT, YY]=funcion(R(ind1,:));
     X=[TT,YY];
     f(ind1,:)=MahalanobisDistance(L,YY);
     [Z, ind]=max(TT);
     z(ind1,:)=YY(ind);
     else
         f(ind1,:)=10;
         z(ind1,:)=0.001;
     end
end
%%%% Encontrar optimo
[glob_op,ind]=min(f);
rB=B(ind,:);
glob_opt=z(ind);
if z(ind,:)>v_opt
    [TT, YY]=funcion(R(ind,:));
    L=[YY];
    v_opt=z(ind,:);
end

%Inicio de iteraciones
while K>k
for ind2=1:Np
pad=5;
al=randperm(Np);
a1=al(1);
a2=al(2);
ri(1,:)=B(ind2,1:150);
r1(1,:)=B(a1,1:150); 
r2(1,:)=B(a2,1:150);
zz=zeros(1,150);
padre=[ri;r1;r2;rB;zz];
hijo1=[];
hijo2=[];
hijo3=[];
hijo4=[];

all=randperm(pad);
rd=rand();
hijo1=[padre(1,1:5+floor(rd*12)),padre(all(1),6+floor(rd*12):6+floor(rd*12)+2),padre(1,7+floor(rd*12)+2:21),padre(1,22:22+floor(rd*17)),padre(all(1),23+floor(rd*17):23+floor(rd*17)+2),padre(1,24+floor(rd*17)+2:42),padre(1,43:43+floor(rd*17)),padre(all(1),44+floor(rd*17):44+floor(rd*17)+2),padre(1,45+floor(rd*17)+2:63),padre(1,64:64+floor(rd*17)),padre(all(1),65+floor(rd*17):65+floor(rd*17)+2),padre(1,66+floor(rd*17)+2:84),padre(1,85:96+floor(rd*6)),padre(all(1),97+floor(rd*6):97+floor(rd*6)+2),padre(1,98+floor(rd*6)+2:105),padre(1,106:106+floor(rd*17)),padre(all(1),107+floor(rd*17):107+floor(rd*17)+2),padre(1,108+floor(rd*17)+2:126),padre(1,127:138+floor(rd*7)),padre(all(1),139+floor(rd*7):139+floor(rd*7)+2),padre(1,140+floor(rd*7)+2:150)];
hijo2=[padre(2,1:5+floor(rd*12)),padre(all(2),6+floor(rd*12):6+floor(rd*12)+2),padre(2,7+floor(rd*12)+2:21),padre(2,22:22+floor(rd*17)),padre(all(2),23+floor(rd*17):23+floor(rd*17)+2),padre(2,24+floor(rd*17)+2:42),padre(2,43:43+floor(rd*17)),padre(all(2),44+floor(rd*17):44+floor(rd*17)+2),padre(2,45+floor(rd*17)+2:63),padre(2,64:64+floor(rd*17)),padre(all(2),65+floor(rd*17):65+floor(rd*17)+2),padre(2,66+floor(rd*17)+2:84),padre(2,85:96+floor(rd*6)),padre(all(2),97+floor(rd*6):97+floor(rd*6)+2),padre(2,98+floor(rd*6)+2:105),padre(2,106:106+floor(rd*17)),padre(all(2),107+floor(rd*17):107+floor(rd*17)+2),padre(2,108+floor(rd*17)+2:126),padre(2,127:138+floor(rd*7)),padre(all(2),139+floor(rd*7):139+floor(rd*7)+2),padre(2,140+floor(rd*7)+2:150)];
hijo3=[padre(3,1:5+floor(rd*12)),padre(all(3),6+floor(rd*12):6+floor(rd*12)+2),padre(3,7+floor(rd*12)+2:21),padre(3,22:22+floor(rd*17)),padre(all(3),23+floor(rd*17):23+floor(rd*17)+2),padre(3,24+floor(rd*17)+2:42),padre(3,43:43+floor(rd*17)),padre(all(3),44+floor(rd*17):44+floor(rd*17)+2),padre(3,45+floor(rd*17)+2:63),padre(3,64:64+floor(rd*17)),padre(all(3),65+floor(rd*17):65+floor(rd*17)+2),padre(3,66+floor(rd*17)+2:84),padre(3,85:96+floor(rd*6)),padre(all(3),97+floor(rd*6):97+floor(rd*6)+2),padre(3,98+floor(rd*6)+2:105),padre(3,106:106+floor(rd*17)),padre(all(3),107+floor(rd*17):107+floor(rd*17)+2),padre(3,108+floor(rd*17)+2:126),padre(3,127:138+floor(rd*7)),padre(all(3),139+floor(rd*7):139+floor(rd*7)+2),padre(3,140+floor(rd*7)+2:150)];
hijo4=[padre(4,1:5+floor(rd*12)),padre(all(4),6+floor(rd*12):6+floor(rd*12)+2),padre(4,7+floor(rd*12)+2:21),padre(4,22:22+floor(rd*17)),padre(all(4),23+floor(rd*17):23+floor(rd*17)+2),padre(4,24+floor(rd*17)+2:42),padre(4,43:43+floor(rd*17)),padre(all(4),44+floor(rd*17):44+floor(rd*17)+2),padre(4,45+floor(rd*17)+2:63),padre(4,64:64+floor(rd*17)),padre(all(4),65+floor(rd*17):65+floor(rd*17)+2),padre(4,66+floor(rd*17)+2:84),padre(4,85:96+floor(rd*6)),padre(all(4),97+floor(rd*6):97+floor(rd*6)+2),padre(4,98+floor(rd*6)+2:105),padre(4,106:106+floor(rd*17)),padre(all(4),107+floor(rd*17):107+floor(rd*17)+2),padre(4,108+floor(rd*17)+2:126),padre(4,127:138+floor(rd*7)),padre(all(4),139+floor(rd*7):139+floor(rd*7)+2),padre(4,140+floor(rd*7)+2:150)];
hijo=[hijo1;hijo2; hijo3;hijo4];
N=[];
ni(1,:)=hijo(1,:);
ni(2,:)=hijo(2,:);
ni(3,:)=hijo(3,:);
ni(4,:)=hijo(4,:);
r11=[1:150];
for ind1=1:4
rand1=rand();
rand2=rand();
if ni(ind1,r11(ceil(rand1*150)))==1
    ni(ind1,(ceil(rand1*150)))=0;
else
    ni(ind1,(ceil(rand1*150)))=1;
end 
if ni(ind1,r11(ceil(rand2*150)))==1
    ni(ind1,(ceil(rand2*150)))=0;
else
    ni(ind1,(ceil(rand2*150)))=1;
end 
end            

N=[hijo;ni];
for ind=1:8
D(ind,:)=[DECOD(N(ind,1:21), info(1,2), info(1,3), info(1,4)), DECOD(N(ind,22:42),info(2,2),info(2,3),info(2,4)), DECOD(N(ind,(43:63)),info(3,2),info(3,3),info(3,4)),DECOD(N(ind,(64:84)),info(4,2),info(4,3),info(4,4)),DECOD(N(ind,(85:105)),info(5,2),info(5,3),info(5,4)),DECOD(N(ind,(106:126)),info(6,2),info(6,3),info(6,4)),DECOD(N(ind,(127:147)),info(7,2),info(7,3),info(7,4)),DECOD(N(ind,(148:149)),info(8,2),info(8,3),info(8,4)),DECOD(N(ind,150),info(9,2),info(9,3),info(9,4))];
end

for ind1=1:8
     Co=funcion2(D(ind1,:));
     if Co<9
     [TT, YY]=funcion(D(ind1,:));
     F(ind1,:)=MahalanobisDistance(L,YY);
     [GG, ind]=max(TT);
     Z(ind1,:)=YY(ind);
     else
         F(ind1,:)=10;
         Z(ind1,:)=0.001;
     end
end
[loc_opt1(ind2), ind]= min(F);
if loc_opt1(ind2)>f(ind2)
    loc_opt1(ind2)=f(ind2);
    B_loc(ind2,:)=B(ind2,:);
    D_loc(ind2,:)=R(ind2,:);
    z_loc(ind2,:)=z(ind2,:);
else
    B_loc(ind2,:)=N(ind,:);
    D_loc(ind2,:)=D(ind,:);
    z_loc(ind2,:)=Z(ind,:);
    
end
end
[glob_o, ind]= min(loc_opt1);
if glob_op>glob_o
    glob_op=glob_o;
    B_glob=B_loc(ind,:);
    D_glob=D_loc(ind,:);
    z_glob=z_loc(ind,:);
   if  glob_opt<z_glob
       glob_opt=z_glob;
   end
if z_glob>v_opt
    [TT, YY]=funcion(D_loc(ind,:));
    L=[YY];
    v_opt=z_loc(ind,:);
    for ind1=1:Np
     Co=funcion2(D_loc(ind1,:));
     if Co<9
     [TT, YY]=funcion(D_loc(ind1,:));
     loc_opt1(ind1)=MahalanobisDistance(L,YY);
     [GG, ind]=max(TT);
     z_loc(ind1,:)=YY(ind);
     else
         loc_opt1(ind1)=10;
         z_loc(ind1,:)=0.001;
     end
    end
end
end
G=sum(-loc_opt1);
G=-loc_opt1./G;
q=[];
for c=1:Np
    q(c,:)=[sum(G(1,1:c))];
end
G=q;
for c=1:Np
padre1(c)=RULETA(G);
end
for ind2=1:Np
pad=5;
al=randperm(Np);
a1=al(1);
a2=al(2);
ri(1,:)=B_loc(padre1(ind2),1:150);
r1(1,:)=B_loc(a1,1:150); 
r2(1,:)=B_loc(a2,1:150);
zz=zeros(1,150);
padre=[ri;r1;r2;rB;zz];
hijo1=[];
hijo2=[];
hijo3=[];
hijo4=[];

all=randperm(pad);
rd=rand();
hijo1=[padre(1,1:5+floor(rd*12)),padre(all(1),6+floor(rd*12):6+floor(rd*12)+2),padre(1,7+floor(rd*12)+2:21),padre(1,22:22+floor(rd*17)),padre(all(1),23+floor(rd*17):23+floor(rd*17)+2),padre(1,24+floor(rd*17)+2:42),padre(1,43:43+floor(rd*17)),padre(all(1),44+floor(rd*17):44+floor(rd*17)+2),padre(1,45+floor(rd*17)+2:63),padre(1,64:64+floor(rd*17)),padre(all(1),65+floor(rd*17):65+floor(rd*17)+2),padre(1,66+floor(rd*17)+2:84),padre(1,85:96+floor(rd*6)),padre(all(1),97+floor(rd*6):97+floor(rd*6)+2),padre(1,98+floor(rd*6)+2:105),padre(1,106:106+floor(rd*17)),padre(all(1),107+floor(rd*17):107+floor(rd*17)+2),padre(1,108+floor(rd*17)+2:126),padre(1,127:138+floor(rd*7)),padre(all(1),139+floor(rd*7):139+floor(rd*7)+2),padre(1,140+floor(rd*7)+2:150)];
hijo2=[padre(2,1:5+floor(rd*12)),padre(all(2),6+floor(rd*12):6+floor(rd*12)+2),padre(2,7+floor(rd*12)+2:21),padre(2,22:22+floor(rd*17)),padre(all(2),23+floor(rd*17):23+floor(rd*17)+2),padre(2,24+floor(rd*17)+2:42),padre(2,43:43+floor(rd*17)),padre(all(2),44+floor(rd*17):44+floor(rd*17)+2),padre(2,45+floor(rd*17)+2:63),padre(2,64:64+floor(rd*17)),padre(all(2),65+floor(rd*17):65+floor(rd*17)+2),padre(2,66+floor(rd*17)+2:84),padre(2,85:96+floor(rd*6)),padre(all(2),97+floor(rd*6):97+floor(rd*6)+2),padre(2,98+floor(rd*6)+2:105),padre(2,106:106+floor(rd*17)),padre(all(2),107+floor(rd*17):107+floor(rd*17)+2),padre(2,108+floor(rd*17)+2:126),padre(2,127:138+floor(rd*7)),padre(all(2),139+floor(rd*7):139+floor(rd*7)+2),padre(2,140+floor(rd*7)+2:150)];
hijo3=[padre(3,1:5+floor(rd*12)),padre(all(3),6+floor(rd*12):6+floor(rd*12)+2),padre(3,7+floor(rd*12)+2:21),padre(3,22:22+floor(rd*17)),padre(all(3),23+floor(rd*17):23+floor(rd*17)+2),padre(3,24+floor(rd*17)+2:42),padre(3,43:43+floor(rd*17)),padre(all(3),44+floor(rd*17):44+floor(rd*17)+2),padre(3,45+floor(rd*17)+2:63),padre(3,64:64+floor(rd*17)),padre(all(3),65+floor(rd*17):65+floor(rd*17)+2),padre(3,66+floor(rd*17)+2:84),padre(3,85:96+floor(rd*6)),padre(all(3),97+floor(rd*6):97+floor(rd*6)+2),padre(3,98+floor(rd*6)+2:105),padre(3,106:106+floor(rd*17)),padre(all(3),107+floor(rd*17):107+floor(rd*17)+2),padre(3,108+floor(rd*17)+2:126),padre(3,127:138+floor(rd*7)),padre(all(3),139+floor(rd*7):139+floor(rd*7)+2),padre(3,140+floor(rd*7)+2:150)];
hijo4=[padre(4,1:5+floor(rd*12)),padre(all(4),6+floor(rd*12):6+floor(rd*12)+2),padre(4,7+floor(rd*12)+2:21),padre(4,22:22+floor(rd*17)),padre(all(4),23+floor(rd*17):23+floor(rd*17)+2),padre(4,24+floor(rd*17)+2:42),padre(4,43:43+floor(rd*17)),padre(all(4),44+floor(rd*17):44+floor(rd*17)+2),padre(4,45+floor(rd*17)+2:63),padre(4,64:64+floor(rd*17)),padre(all(4),65+floor(rd*17):65+floor(rd*17)+2),padre(4,66+floor(rd*17)+2:84),padre(4,85:96+floor(rd*6)),padre(all(4),97+floor(rd*6):97+floor(rd*6)+2),padre(4,98+floor(rd*6)+2:105),padre(4,106:106+floor(rd*17)),padre(all(4),107+floor(rd*17):107+floor(rd*17)+2),padre(4,108+floor(rd*17)+2:126),padre(4,127:138+floor(rd*7)),padre(all(4),139+floor(rd*7):139+floor(rd*7)+2),padre(4,140+floor(rd*7)+2:150)];
hijo=[hijo1;hijo2; hijo3;hijo4];
N=[];
ni(1,:)=hijo(1,:);
ni(2,:)=hijo(2,:);
ni(3,:)=hijo(3,:);
ni(4,:)=hijo(4,:);
r11=[1:150];
for ind1=1:4
rand1=rand();
rand2=rand();
if ni(ind1,r11(ceil(rand1*150)))==1
    ni(ind1,(ceil(rand1*150)))=0;
else
    ni(ind1,(ceil(rand1*150)))=1;
end 
if ni(ind1,r11(ceil(rand2*150)))==1
    ni(ind1,(ceil(rand2*150)))=0;
else
    ni(ind1,(ceil(rand2*150)))=1;
end 
end            

N=[hijo;ni];
for ind=1:8
D(ind,:)=[DECOD(N(ind,1:21), info(1,2), info(1,3), info(1,4)), DECOD(N(ind,22:42),info(2,2),info(2,3),info(2,4)), DECOD(N(ind,(43:63)),info(3,2),info(3,3),info(3,4)),DECOD(N(ind,(64:84)),info(4,2),info(4,3),info(4,4)),DECOD(N(ind,(85:105)),info(5,2),info(5,3),info(5,4)),DECOD(N(ind,(106:126)),info(6,2),info(6,3),info(6,4)),DECOD(N(ind,(127:147)),info(7,2),info(7,3),info(7,4)),DECOD(N(ind,(148:149)),info(8,2),info(8,3),info(8,4)),DECOD(N(ind,150),info(9,2),info(9,3),info(9,4))];
end
Co=funcion2(D(ind1,:));
     if Co<9
     [TT, YY]=funcion(D(ind1,:));
     F(ind1,:)=MahalanobisDistance(L,YY);
     [GG, ind]=max(TT);
     Z(ind1,:)=YY(ind);
     else
         F(ind1,1)=10;
         Z(ind1,1)=0.001;
     end
[loc(ind2), ind]= min(F);
if loc(ind2)>loc_opt1(ind2)
    loc(ind2)=loc_opt1(ind2);
    Bloc(ind2,:)=B_loc(ind2,:);
    Dloc(ind2,:)=D_loc(ind2,:);
    zloc(ind2,1)=z_loc(ind2,1);
else
    Bloc(ind2,:)=N(ind,:);
    Dloc(ind2,:)=D(ind,:);
    zloc(ind2,:)=Z(ind,:);
    B_loc(ind2,:)=N(ind,:);
    end
end
opti=[loc_opt1,loc];
B_opt=[Bloc;B_loc];
D_opt=[Dloc;D_loc];
z_opt=[z_loc;zloc];
[glob_opt1, ind]= max(z_opt);
if glob_opt<glob_opt1
    glob_opt=glob_opt1;
    B_glob=B_opt(ind,:);
    D_glob=D_opt(ind,:);
    g_opt=opti(ind);
if glob_opt>v_opt
    [TT, YY]=funcion(D_opt(ind,:));
    L=[YY];
    v_opt=glob_opt;
    for ind1=1:Np
     Co=funcion2(D_loc(ind1,:));
     if Co<9
     [TT, YY]=funcion(D_loc(ind1,:));
     loc_opt1(ind1)=MahalanobisDistance(L,YY);
     [GG, ind]=max(TT);
     z_loc(ind1,1)=YY(ind);
     else
         loc_opt1(1,ind1)=10;
         z_loc(ind1,1)=YY(ind);
     end
    end
end
end
mejor(k)=glob_opt;
disp(sprintf('generacion = %d, mejor fenotipo = %.6f', k, glob_opt))
k=k+1;
end
EXPERIMENTOS(ind3,:)=[D_glob,glob_opt];
SALIDA=EXPERIMENTOS(ind3,:);
save('EXPERIMENTO_1a.txt','SALIDA','-ASCII','-append');
end
k=k;
% funcion1(D_glob);
% figure
% plot(mejor);
% title('Convergencia')
% ylabel('angulo')
% xlabel('iteracion')