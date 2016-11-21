function learnEfficiency

data = csvread('efficiencies.csv');

horsepower = data(:,1);
mass = data(:,3);
%targetA = data(:,6);
%targetB = data(:,7);
%targetC = data(:,8);

efficiencies = data(:,11);

N = length(mass);
%{
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

test_hp = horsepower(testing); train_hp = horsepower(training);
test_mass = mass(testing); train_mass = mass(training);
test_eff = efficiencies(testing); train_eff = efficiencies(training);
%}

% let's try KNN as a concept
K=3;
prediction = zeros(N,1);

for i = 1:N
    distances = zeros(N,2);
    
    for j = 1:N
        diff = [mass(i);horsepower(i)]-[mass(j);horsepower(j)];
        distances(j,:) = [norm(diff,2),j];
    end

    distances = sortrows(distances);
    % find k closest points
    index = distances(1:K,2);

    f = 0;

    for j = 1:K
        f = f + efficiencies(index(j));
    end

    prediction(i) = f/K;
end

figure(1)
plot([1:N],[efficiencies,prediction])

av_error = sum(abs(efficiencies-prediction))/N
end


%{
% for each testing point
for i = 1:length(testing)
    distances = zeros(length(training),2);
    
    %compute distance to every training point
    for j = 1:length(training)
        distances(j,:) = [norm([test_hp(i);test_mass(i)]-...
            [train_hp(j);train_mass(j)],2),j];
    end
%}    
    

