% This is a function to specify growth and secretion parameters given
% specific regulatory mechanisms.
% Here we only specify parameters that vary with regulation mechanism.
% common growth parameters (e.g. gmax will be specified in the main script
% file)

function [secparams,gparams,decayparams] = SpecifyRuleParams(ruletype,effparams)

if strcmp(ruletype,'1ginhib')
    % Here we just have 1 growth inhibitor that is secreted by all cells at
    % a constant rate (no regulation)
    
    % extract parameters from effective parameters
    mu = effparams.mueff;
    
    % secretion rules
    secparams.mumaxvec = mu;
    secparams.ifmumaxreg = false;
    secparams.ifregsec = false;
    
    % decay parameters
    decayparams.decayRates = 1;
    
    % growth parameters
    gparams.nvec = -100;
    gparams.ifreggthres = false;

elseif strcmp(ruletype,'1ginhib_1gthresreg')
    % Here we have 1 growth inhibitor that is secreted by all cells at a
    % constant rate, and the second chemical is a growth-threshold
    % regulator that is secreted only when concentration of the growth
    % inhibitor is below or above some threshold. The maximum secretion
    % rate of the growth threshold regulator can also be regulated by c1.
    
    % extract parameters from effective parameters
    mu1eff = effparams.mu1eff;
    mu2eff = effparams.mu2eff;
    k2overk1 = effparams.k2overk1;
    mu2max_regcoeff = effparams.mu2max_regcoeff;
    
    secthres_regsign = effparams.secthres_regsign; % must be either 1 or -1 or 0
    if secthres_regsign ~= 0
        secthres_eff = effparams.secthres_eff;
    end
    gthres_regsign = effparams.gthres_regsign; % must be either 1 or -1
    
    % secretion rules
    secparams.mumaxvec = [mu1eff;mu2eff];
    if mu2max_regcoeff == 0
        secparams.ifmumaxreg = false;
    else
        secparams.ifmumaxreg = true;
        secparams.muregmat = [0, 0; mu2max_regcoeff, 0];
        secparams.muregpower = 1;
    end
    if secthres_regsign == 0
        secparams.ifregsec = false;
    elseif secthres_regsign == -1 || secthres_regsign == 1 
        secparams.ifregsec = true;
        secparams.nijmat = [0, 0; secthres_regsign*100, 0];
        secparams.Kijmat = [0, 0; secthres_eff, 0];
        secparams.ifregsecthres = false;
    else
        error('secthres_regsign must be either 0, 1 or -1!');
    end
    
    % decay parameters
    decayparams.decayRates = [1;k2overk1];
    
    % growth parameters
    gparams.nvec = [-100,0];
    gparams.ifreggthres = true;
    gparams.kMat = [0, gthres_regsign; 0, 0];
    gparams.kpower = 1;
    
elseif strcmp(ruletype,'2ginhib')
    % Here we have 2 growth inhibitors, the first uniformly secreted while
    % the second is only secreted when the concentration of the first is
    % below some threshold. The maximum secretion rate of the second growth
    % regulator can also be a function of c1.

    % extract parameters from effective parameters
    mu1eff = effparams.mu1eff;
    mu2eff = effparams.mu2eff;
    k2overk1 = effparams.k2overk1;
    secthres_eff = effparams.secthres_eff;
    mu2max_regcoeff = effparams.mu2max_regcoeff;
    
    % secretion rules
    secparams.mumaxvec = [mu1eff;mu2eff];
    if mu2max_regcoeff == 0
        secparams.ifmumaxreg = false;
    else
        secparams.ifmumaxreg = true;
        secparams.muregmat = [0, 0; mu2max_regcoeff, 0];
        secparams.muregpower = 1;
    end
    secparams.ifregsec = true;
    secparams.nijmat = [0, 0; -100, 0];
    secparams.Kijmat = [0, 0; secthres_eff, 0];
    secparams.ifregsecthres = false;
    
    % decay parameters
    decayparams.decayRates = [1;k2overk1];
    
    % growth parameters
    gparams.nvec = [-100,-100];
    gparams.ifreggthres = false;
    
end

end











