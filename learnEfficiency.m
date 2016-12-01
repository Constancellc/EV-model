function learnEfficiency

y = csvread('efficiencies.csv');

data = csvread('../EPA-Code/pureEVdata.csv',0,3);

%{
 1 |   horsepower    |
 2 | # driven wheels |    -
 3 |   kerb weight   |   lbs
 4 |    axle ratio   |    -
 5 |    n/v ratio    |    -
 6 |  Target coeff A |   lbf
 7 |  Target coeff B |  lbf/mph
 8 |  Target coeff C | lbf/mph^2
 9 |   Highways fe   |   MPGe
10 |     UDDS fe     |   MPGe

%}
%{
figure(1)
for i = 1:8
    subplot(2,4,i)
    plot(data(:,i),efficiencies,'x')
end
%}

% let's assume vehicles are defined by horsepower, weight, axle ratio

X = [data(:,1),data(:,3),data(:,4)];

N = length(X);

% {
testing = [];
training = [];

for i = 1:N
    n = rand;
    
    if n > 0.7
        testing = [testing;i];
    else
        training = [training;i];
    end
end

X_training = X(training,:); X_testing = X(testing,:)
y_training = y(training,:); y_testing = y(testing,:);
%}


% {
% let's try KNN as a concept
K=2;

M = length(testing);
prediction = zeros(M,1);

for i = 1:M
    distance = zeros(N-M,2);
    
    for j = 1:(N-M)
        gap = X_training(j,:)-X_testing(i,:)
        distance(j,:) = [norm(gap,2),j];
    end

    distance = sortrows(distance)
    
    % find k closest points
    index = distance(1:K,2)

    f = 0;

    for j = 1:K
        f = f + y_training(index(j));
    end

    prediction(i) = f/K;
end

figure(1)
plot([1:M],[prediction,y_testing],'x')