function C=Covariance(X)
[n,k]=size(X);
Xc=X-repmat(mean(X),n,1);
C=Xc'*Xc/n;
end