function  [X_robust,Pend,iter] = robust(A,L,P,k1,k2,numJ,e)
format long
[n,t] = size(A);
X1 = (A'*P*A)\(A'*P*L);               
X2 = 99999*ones(t,1);                 
V = L - A*X1;                         
sigma0 = sqrt((V'*P*V)/(n-t));         
% PV = inv(P)- A*inv(A'*P*A)*A';  
PV = diag(1./diag(P)) - A*((A'*P*A)\A');
PVV = sigma0*sqrt(diag(PV));          
iter = 0;
while max(abs(X2-X1))>e/1000 && iter <15
      
      iter = iter + 1;
      X2 = X1;
      V = A*X2 - L;                     
      abVV = abs(V./PVV);               
      tao =  down_t(abVV,k1,k2,numJ);   
      down_P = diag(tao.*diag(P));      
      X1 = (A'*down_P*A)\(A'*down_P*L); 
end
X_robust = X1;       
Pend = down_P;        