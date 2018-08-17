% Alex Williams Alternate Minimization algorithm for PCA
K=2; % number of components
numIterations=1000; % number of iterations to perform
%data = randn(100,K) * randn(K, 101);

data=myOR';
[M, N] = size(data);
U = randn(M,K);
loss=zeros(1,numIterations);

for iteration = 1:numIterations
   Vt = U \ data; % update V (fixed U)
   U = data / Vt; % update U (fixed V)
   loss(iteration) = norm(data - U*Vt,'fro');
   
   if mod(iteration,10)==0
      display(['iteration #' num2str(iteration)]) 
   end
end