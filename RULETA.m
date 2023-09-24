function seleccionado = RULETA(E)
r=rand();
seleccionado = -1;
for ind=1:length(E)
    if(E(ind)>r)
        seleccionado = ind;
        break;
    end
end
end
