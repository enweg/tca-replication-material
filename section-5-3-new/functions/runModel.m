function results = runModel(modelPath, modelName, maxHorizon)
    global M_ options_ oo_
    
    results = struct(); 
    originalDir = pwd; 
    cd(modelPath); dynare(modelName); cd(originalDir);

    results.M_ = M_; 
    results.options_ = options_; 
    results.oo_ = oo_; 

    % --- TOTAL EFFECTS
    model = DSGE(M_, options_, oo_);
    results.irfsVarma = model.IRF(maxHorizon).irfs; 
    results.irfsSS = model.IRFStateSpace(maxHorizon).irfs; 
end
