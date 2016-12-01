N = 100;

errors = zeros(4,1);

for i = 1:N
    
    rmse = modelTesting;
    
    errors = errors + rmse;
end

errors = errors/N