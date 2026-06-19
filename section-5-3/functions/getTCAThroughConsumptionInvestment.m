% Channels are defined as going from consumption directly to output-inflation 
% and not through investment. Similarly for investment channel. 
% More precisely, let M_{i, j} be the set of paths that propagate within 
% interest rates and consumption up to period i, leave consumption in period i 
% and directly go into output-inflation system in period j. The consumption 
% channel in period h then consists of the union of M_{i,j} for i<=j<=h. 


function results = getTCAThroughConsumptionInvestment(runOptions, results)
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
            MijConsumption = makeMijConsumption(i, j, order); 
            MijInvestment= makeMijInvestment(i, j, order); 

            eConsumption = model.transmission('eps_nu', MijConsumption, order, maxHorizon); 
            eInvestment = model.transmission('eps_nu', MijInvestment, order, maxHorizon); 

            effectsConsumption(:, :, (j + 1):end) = effectsConsumption(:, :, (j + 1):end) + eConsumption(:, :, (j + 1):end); 
            effectsInvestment(:, :, (j + 1):end) = effectsInvestment(:, :, (j + 1):end) + eInvestment(:, :, (j + 1):end); 
        end
    end

    if ~isfield(results, 'TCAThroughConsumptionInvestment')
        results.TCAThroughConsumptionInvestment = struct(); 
        results.TCAThroughConsumptionInvestment.orderCBeforeITrue = struct(); 
        results.TCAThroughConsumptionInvestment.orderCBeforeIFalse = struct(); 
    end

    if runOptions.orderCBeforeI
        results.TCAThroughConsumptionInvestment.orderCBeforeITrue.effectsConsumption = effectsConsumption; 
        results.TCAThroughConsumptionInvestment.orderCBeforeITrue.effectsInvestment = effectsInvestment; 
    else
        results.TCAThroughConsumptionInvestment.orderCBeforeIFalse.effectsConsumption = effectsConsumption; 
        results.TCAThroughConsumptionInvestment.orderCBeforeIFalse.effectsInvestment = effectsInvestment; 
    end

end

function Q = makeMijConsumption(i, j, order)
    assert(i <= j); 

    % original order is c, I, pi
    Q = through(1, i, order) & through(3, j, order); 
    if i > 0
        Q = Q & notThrough(2:3, 0:(i - 1), order); 
    end
    if j == i
        Q = Q & notThrough(2, i, order); 
    elseif j == i + 1
        Q = Q & notThrough(2:3, i, order) & notThrough(1:2, j, order); 
    else
        Q = Q & notThrough(2:3, i, order) & notThrough(1:3, (i + 1):(j - 1), order) & notThrough(1:2, j, order);
    end
end

function Q =  makeMijInvestment(i, j, order)
    assert(i <= j); 

    % original order is c, I, pi
    Q = through(2, i, order) & through(3, j, order); 
    if i > 0
        Q = Q & notThrough([1, 3], 0:(i - 1), order); 
    end
    if j == i
        Q = Q & notThrough(1, i, order); 
    elseif j == i + 1
        Q = Q & notThrough([1, 3], i, order) & notThrough(1:2, j, order);
    else
        Q = Q & notThrough([1, 3], i, order) & notThrough(1:3, (i + 1):(j - 1), order) & notThrough(1:2, j, order); 
    end
end


