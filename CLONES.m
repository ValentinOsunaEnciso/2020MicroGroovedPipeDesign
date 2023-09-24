
% Reproduction
function [T,pcs] = CLONES(n,Pc,N,ind,P)
% n		-> number of clones
% fat	-> multiplying factor
% ind	-> best individuals
% T		-> temporary population
% pcs	-> final position of each clone
T=[]; 
for i1=1:n
      cs(i1) = round(Pc*N);
      pcs(i1) = sum(cs);
      T = [T; ones(cs(i1),1) * P(ind(end-i1+1),:)];
end

