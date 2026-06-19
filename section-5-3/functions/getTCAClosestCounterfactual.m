% The closest we can get to the counterfactual setup in which consumption and 
% investment are not moving at all is to take a channel which never goes through 
% consumption and investment respectively. 

function results = getTCAClosestCounterfactual(runOptions, results)
    if ~isfield(runOptions, 'orderCBeforeI')
        runOptions.orderCBeforeI = true; 
    end
    if ~isfield(runOptions, 'maxHorizon')
        error("maxHorizon option must be specified");
    end
    maxHorizon = runOptions.maxHorizon; 

    if runOptions.orderCBeforeI
        order = {'c', 'I', 'pi'};
    else
        order = {'I', 'c', 'pi'};
    end

    model = DSGE(results.M_, results.options_, results.oo_);
    effectNotConsumption = zeros(3, 1, maxHorizon + 1); 
    effectNotInvestment = zeros(3, 1, maxHorizon + 1);

    for h = 0:maxHorizon
        QConsumption = model.notThrough('c', 0:h, order);  
        QInvestment = model.notThrough('I', 0:h, order); 

        eConsumption = model.transmission('eps_nu', QConsumption, order, maxHorizon); 
        eInvestment = model.transmission('eps_nu', QInvestment, order, maxHorizon); 

        effectNotConsumption(:, :, h + 1) = eConsumption(:, :, h + 1); 
        effectNotInvestment(:, :, h + 1) = eInvestment(:, :, h + 1); 
    end

    if ~isfield(results, 'TCAClosestCounterfactual')
        results.TCAClosestCounterfactual = struct(); 
        results.TCAClosestCounterfactual.orderCBeforeITrue = struct(); 
        results.TCAClosestCounterfactual.orderCBeforeIFalse = struct();
    end

    if runOptions.orderCBeforeI
        results.TCAClosestCounterfactual.orderCBeforeITrue.effectNotConsumption = effectNotConsumption; 
        results.TCAClosestCounterfactual.orderCBeforeITrue.effectNotInvestment = effectNotInvestment; 
    else
        results.TCAClosestCounterfactual.orderCBeforeIFalse.effectNotConsumption = effectNotConsumption; 
        results.TCAClosestCounterfactual.orderCBeforeIFalse.effectNotInvestment = effectNotInvestment; 
    end

end
