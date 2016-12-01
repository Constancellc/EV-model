function effPred = KNN(X_training,X_testing,y)

K = 3;

N = length(X_training); M = length(X_testing);

effPred = zeros(M,1);

for i = 1:M
    distances = zeros(N,2);
    
    for j = 1:N
        distances(j,:) = [norm(X_training(j,:)-X_testing(i,:),2),j];
    end
    
    distances = sortrows(distances);
    
    f = 0;
    
    for j = 1:K
        f = f + y(distances(j,2));
    end
    
    effPred(i) = f/K;
    
end

