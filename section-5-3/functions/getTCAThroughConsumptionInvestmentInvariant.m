% Computes consumption and investment channel. 
% Effects should be invariant to re-ordering
% TODO: prove this. 

function results = getTCAThroughConsumptionInvestmentSensitive(runOptions, results)
    if ~isfield(runOptions, 'orderCBeforeI')
        runOptions.orderCBeforeI = true; 
    end
    if ~isfield(runOptions, 'maxHorizon')
        error("maxHorizon option must be specified");
    end
    maxHorizon = runOptions.maxHorizon; 

    model = DSGE(results.M_, results.options_, results.oo_); 

    % Check if observation order equals R, c, I, Y, pi
    if ~isequal(DSGE.defineOrder_(model.getVariableNames(), results.options_)', 1:3)
        error("Observation order should be c, I, pi");
    end

    if runOptions.orderCBeforeI
        orderVars = {'c', 'I', 'pi'};
        order =  model.defineOrder(orderVars); 
        if any(order ~= (1:3)')
            error("Order not supported")
        end
    else
        orderVars = {'I', 'c', 'pi'}
        order = model.defineOrder(orderVars); 
        if any(order ~= [2; 1; 3])
            error("Order not supported")
        end
    end



    effectsConsumption = zeros(3, 1, maxHorizon + 1); 
    effectsInvestment = zeros(3, 1, maxHorizon + 1);
    for j = 0:maxHorizon
        for i = 0:j
            MijConsumption = makeMijConsumption(i, j, order, runOptions.orderCBeforeI); 
            MijInvestment= makeMijInvestment(i, j, order, runOptions.orderCBeforeI); 

            eConsumption = model.transmission('eps_nu', MijConsumption, order, maxHorizon); 
            eInvestment = model.transmission('eps_nu', MijInvestment, order, maxHorizon); 

            effectsConsumption(:, :, (j + 1):end) = effectsConsumption(:, :, (j + 1):end) + eConsumption(:, :, (j + 1):end); 
            effectsInvestment(:, :, (j + 1):end) = effectsInvestment(:, :, (j + 1):end) + eInvestment(:, :, (j + 1):end); 
        end
    end

    if ~isfield(results, 'TCAThroughConsumptionInvestmentInvariant')
        results.TCAThroughConsumptionInvestmentInvariant = struct(); 
        results.TCAThroughConsumptionInvestmentInvariant.orderCBeforeITrue = struct(); 
        results.TCAThroughConsumptionInvestmentInvariant.orderCBeforeIFalse = struct(); 
    end

    if runOptions.orderCBeforeI
        results.TCAThroughConsumptionInvestmentInvariant.orderCBeforeITrue.effectsConsumption = effectsConsumption; 
        results.TCAThroughConsumptionInvestmentInvariant.orderCBeforeITrue.effectsInvestment = effectsInvestment; 
    else
        results.TCAThroughConsumptionInvestmentInvariant.orderCBeforeIFalse.effectsConsumption = effectsConsumption; 
        results.TCAThroughConsumptionInvestmentInvariant.orderCBeforeIFalse.effectsInvestment = effectsInvestment; 
    end

end

function Q = makeMijConsumption(i, j, order, orderCBeforeI)
    assert(i <= j);

    if orderCBeforeI
        Q = makeMijConsumptionCBeforeI(i, j, order);
    else
        Q = makeMijConsumptionIBeforeC(i, j, order);
    end
end
function Q = makeMijInvestment(i, j, order, orderCBeforeI)
    assert(i <= j);

    if orderCBeforeI
        Q = makeMijInvestmentCBeforeI(i, j, order);
    else
        Q = makeMijInvestmentIBeforeC(i, j, order);
    end
end

function Q = makeMijConsumptionCBeforeI(i, j, order)
    assert(i <= j);

    Q = through(1, i, order) & through(3, j, order); 
    if i > 0
        Q = Q & notThrough(3, 0:(i - 1), order);
    end
    if j == i 
        Q = Q & notThrough(2, j, order);
    elseif j == i + 1
        Q = Q & notThrough(2:3, i, order) & notThrough(1:2, j, order); 
    else
        Q = Q & notThrough(2:3, i, order) & notThrough(1:3, (i + 1):(j - 1), order) & notThrough(1:2, j, order);
    end
end
function Q = makeMijConsumptionIBeforeC(i, j, order)
    assert(i <= j);

    Q = through(1, i, order) & through(3, j, order); 
    if i > 0
        Q = Q & notThrough(3, 0:(i - 1), order); 
    end
    if j == i
        % nothing to do
    elseif j == i + 1
        Q = Q & notThrough(3, i, order) & notThrough(1:2, j, order); 
    else 
        Q = Q & notThrough(3, i, order) & notThrough(1:3, (i + 1):(j - 1), order) & notThrough(1:2, j, order);
    end
end

function Q = makeMijInvestmentCBeforeI(i, j, order)
    assert(i <= j);

    Q = through(2, i, order) & through(3, j, order); 
    if i > 0
        Q = Q & notThrough(3, 0:(i - 1), order);
    end
    if j == i
        % Q = Q & notThrough(1, j, order); 
    elseif j == i + 1
        Q = Q & notThrough(3, i, order) & notThrough(1:2, j, order);
    else
        Q = Q & notThrough(3, i, order) & notThrough(1:3, (i + 1):(j - 1), order) & notThrough(1:2, j, order);
    end
end
function Q = makeMijInvestmentIBeforeC(i, j, order)
    assert(i <= j);

    Q = through(2, i, order) & through(3, j, order); 
    if i > 0
        Q = Q & notThrough(3, 0:(i - 1), order);
    end
    if j == i
        Q = Q & notThrough(1, j, order); 
    elseif j == i + 1
        Q = Q & notThrough([1, 3], i, order) & notThrough(1:2, j, order);
    else
        Q = Q & notThrough([1, 3], i, order) & notThrough(1:3, (i + 1):(j - 1), order) & notThrough(1:2, j, order);
    end
end
