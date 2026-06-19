function results = getTotalEffects(runOptions, results)

    if ~isfield(runOptions, 'maxHorizon')
        error('maxHorizon option must be specified')
    end
    maxHorizon = runOptions.maxHorizon; 

    global M_ options_ oo_;
    if nargin < 2
        results = struct();
    end
    
    model = DSGE(M_, options_, oo_);
    results.irfsVarma = model.IRF(maxHorizon).irfs; 
    results.irfsSS = model.IRFStateSpace(maxHorizon).irfs; 
end
