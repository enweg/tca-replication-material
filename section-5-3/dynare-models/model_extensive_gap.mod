%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL NAME: GALI STYLE DSGE WITH INVESTMENT ADJUSTMENT COSTS
% AUTHOR: Enrico Wegner (e.wegner@maastrichtuniversity.nl)
% DATE: 2025/03/06
% DESCRIPTION: Largely uses the setup of Gali 2015 Chapter 3, but includes 
%              investment adjustment costs as in, e.g., Christiano etal. 2005.
%              This version of the model includes the extensive set of
%              equations; i.e. equations have not yet been reduced to a smaller
%              set. Therefore, this file can be used to test the models with 
%              the reduced number of equations. This file also includes the 
%              the natural economy block and computes gaps. 
% REFERENCES: 
% - Christiano, L. J., Eichenbaum, M., & Evans, C. L. (2005). 
%   Nominal Rigidities and the Dynamic Effects of a Shock to Monetary Policy. 
%   Journal of Political Economy, 113(1), 1–45. https://doi.org/10.1086/426038
% - Justiniano, A., Primiceri, G. E., & Tambalotti, A. (2011). 
%   Investment shocks and the relative price of investment. 
%   Review of Economic Dynamics, 14(1), 102–121. 
%   https://doi.org/10.1016/j.red.2010.08.004
% - Galí, J. (2015). Monetary policy, inflation, and the business cycle: 
%   An introduction to the new Keynesian framework and its applications 
%   (Second edition). Princeton University Press.
% - Gertler, M., & Karadi, P. (2011). A model of unconventional monetary policy. 
%   Journal of Monetary Economics, 58(1), 17–34. 
%   https://doi.org/10.1016/j.jmoneco.2010.10.004
% NOTES
% - We denote with K_t capital at the beginning of period t. This is different
%   from the default Dynare implementation. 'predetermined_variables' fixes this.
% - Possibility of interest rate smoothing was implemented to achieve more 
%   plausible IRFs. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/*------------------------------------------------------------------------------
ENDOGENOUS VARIABLES
------------------------------------------------------------------------------*/
var
    c         $\hat{c}$          (long_name='Consumption (sticky-price; deviation)')
    c_n       $\hat{c}^n$        (long_name='Consumption (natural economy; deviation)')
    c_g       $\hat{c}^g$        (long_name='Consumption gap')

    I         $\hat{I}$          (long_name='Investment (sticky-price; deviation)')
    I_n       $\hat{I}^n$        (long_name='Investment (natural economy; deviation)')
    I_g       $\hat{I}^g$        (long_name='Investment gap')

    K         $\hat{K}$          (long_name='Capital (sticky-price; deviation)')
    Keff
    K_n       $\hat{K}^n$        (long_name='Capital (natural economy; deviation')
    Keff_n
    K_g       $\hat{K}^g$        (long_name='Capital gap')

    Y         $\hat{Y}$          (long_name='Output (sticky-price; deviation)')
    Y_n       $\hat{Y}^n$        (long_name='Output (natural economy; deviation)')
    Y_g       $\hat{Y}^g$        (long_name='Output gap')

    R         $\hat{R}$          (long_name='Gross nominal interest on bonds (sticky-price; deviation)')
    r_n       $\hat{r}^n$        (long_name='Real natural rate of interest (deviation)') 
    pi        $\hat{\pi}$        (long_name='Gross inflation (sticky-price; deviation)')

    N         $\hat{N}$          (long_name='Hours worked (sticky-price; deviation)')
    N_n       $\hat{N}^n$        (long_name='Hours worked (natural economy; deviation)')
    N_g       $\hat{N}^g$        (long_name='Employment gap')

    Wr        $\hat{W}^R$        (long_name='Real wage (sticky-price; deviation)')
    Wr_n      $(\hat{W}^R)^n$    (long_name='Real wage (natural economy; deviation)')
    Wr_g      $(\hat{W}^R)^g$    (long_name='Real wage gap')
    
    q         $\hat{q}$          (long_name='Tobins q (sticky-price; deviation)')
    q_n       $\hat{q}^n$        (long_name='Tobins q (natural economy; deviation)')
    q_g       $\hat{q}^g$        (long_name='Tobins q gap')

    rk        $\hat{r}^k$        (long_name='Real rental rate of capital (sticky-price; deviation)')
    rk_n      $(\hat{r}^k)^n$    (long_name='Real rental rate of capital (natural economy; deviation)')
    rk_g      $(\hat{r}^k)^g$    (long_name='Real rental rate of capital gap')

    mc        $\lambda^f$        (long_name='Real marginal cost (sticky-price; deviation)')

    z         $\hat{Z}$          (long_name='Preference shock')
    // pI        $\hat{p}^I$        (long_name='Investment-specific technology (IST) shock')
    // A         $\hat{A}$          (long_name='Technology shock')
    Omega     $\hat{\Omega}$     (long_name='Capital quality shock')
    nu        $\nu$              (long_name='Monetary policy shock')
;

/*------------------------------------------------------------------------------
EXOGENOUS SHOCKS / INNOVATIONS
------------------------------------------------------------------------------*/
varexo
    eps_z     $\varepsilon^Z$        (long_name='Preference shock innovation')
    // eps_pI    $\varepsilon^{P^I}$    (long_name='IST shock innovation')
    // eps_A     $\varepsilon^{A}$      (long_name='Technology shock innovation')
    eps_Omega $\varepsilon^{\Omega}$ (long_name='Capital quality shock innovation')
    eps_nu    $\varepsilon^\nu$      (long_name='Monetary policy shock innovation')
;

/*------------------------------------------------------------------------------
PARAMETERS
------------------------------------------------------------------------------*/
parameters
    beta          $\beta$       (long_name='Discount factor')
    delta         $\delta$      (long_name='Depreciation rate')
    sigma         $\sigma$      (long_name='Coefficient of relative risk aversion')
    varphi        $\varphi$     (long_name='Inverse Frisch elasticity of labour supply')
    rk_bar        $r^k$         (long_name='Steady-state real capital rental rate')
    s_c           $s_C$         (long_name='Steady-state consumption to output (sticky-price)')
    s_I           $s_I$         (long_name='Steady-state investment to output (sticky-price)')
    phi_I         $\phi_I$      (long_name='Investment adjustment cost curvature')
    alpha         $\alpha$      (long_name='Capital share in production')
    theta         $\theta$      (long_name='Price non-adjustment probability')
    epsilon       $\epsilon$    (long_name='Elasticity of substitution across differentiated goods')

    Phi_pi        $\Phi_\pi$    (long_name='Taylor rule coefficient on inflation')
    Phi_y         $\Phi_y$      (long_name='Taylor rule coefficient on output (gap)')

    rho_z         $\rho_Z$      (long_name='Persistence of preference shock')
    // rho_pI        $\rho_{P^I}$  (long_name='Persistence of IST shock')
    // rho_A         $\rho_A$      (long_name='Persistence of technology shock')
    rho_Omega     $\rho_\Omega$ (long_name='Persistence of capital quality shock')
    rho_nu        $\rho_\nu$    (long_name='Persistence of monetary policy shock')
    rho_R         $\rho_R$      (long_name='Interest rate smoothing / persistence')
;


/*------------------------------------------------------------------------------
CALIBRATION
------------------------------------------------------------------------------*/

// Deep structural
beta = 1.03^(-0.25);          // Christiano etal. (2005) p. 15
alpha = 0.36;                 // Christiano etal. (2005) p. 15
delta = 0.025;                // Christiano etal. (2005) p. 15
sigma = 1;                    // Used in both Gali (2015) and Christiano etal. (2005) [They use log utility]
varphi = 5;                   // Gali (2015)
epsilon = 9;                  // Gali (2015)
theta = 3 / 4;                // Gali (2015)
phi_I = 2.48;                 // Christiano etal. (2005) Table 2 [Called kappa there]

// Monetary policy
Phi_pi = 1.5;                 // Gali (2015) p. 68
Phi_y = 0.5 / 4;              // Gali (2015) p. 68
rho_R = 0;                    // Baseline has no interest rate smoothing

// Shocks persistence
rho_z = 0.5;                  // Gali (2015) Section 3.4.1.2
// rho_A = 0.9;                  // Gali (2015) Section 3.4.1.3
// rho_pI = 0.19;                // Justiano etal. (2009)
rho_Omega = 0.66;             // Gertler and Karadi (2011) p. 27 (p. 11 in the pdf)
// rho_nu = 0.5;                 // Gali(2015) p. 68
rho_nu = 0.3;                 // Less than Gali to ensure positive R response

steady_state_model;
    // Derived 
    rk_bar = (1 - beta * (1 - delta)) / beta;
    s_I = alpha * (epsilon - 1) / epsilon * delta * (1 / rk_bar);
    s_c = 1 - s_I; 
end;

/*------------------------------------------------------------------------------
MODEL DEFINITION
- The model is already linearised by hand. 
------------------------------------------------------------------------------*/
// Capital is predetermined in our notation
predetermined_variables K K_n;

model(linear);
    Keff = Omega + K; 
    Keff_n = Omega + K_n;
    /////////////////////////////////////////////////////////////////////////// 
    // --- SHOCK PROCESSES
    /////////////////////////////////////////////////////////////////////////// 
    z = rho_z * z(-1) + eps_z;
    // pI = rho_pI * pI(-1) + eps_pI;
    // A = rho_A * A(-1) + eps_A;
    Omega = rho_Omega * Omega(-1) + eps_Omega;

    /////////////////////////////////////////////////////////////////////////// 
    // --- HOUSEHOLD BLOCK 
    // Sticky equations start with LH, natural economy equations start with NE
    /////////////////////////////////////////////////////////////////////////// 
    // LH1
    c = c(+1) - 1 / sigma * (R - pi(+1)) - 1 / sigma * (rho_z - 1) * z;
    // NE11
    c_n = c_n(+1) - 1 / sigma * r_n - 1 / sigma * (rho_z - 1) * z;

    // LH2
    varphi * N = Wr - sigma * c;
    // NE10
    varphi * N_n = Wr_n - sigma * c_n;

    // LH3
    q = beta * rk_bar * rk(+1) + beta * (1 - delta) * q(+1) - (R - pi(+1));
    // NE9
    q_n = beta * rk_bar * rk_n(+1) + beta * (1 - delta) * q_n(+1) - r_n;

    // LH4
    // q - pI = phi_I * ((I - I(-1)) - beta * (I(+1) - I));
    q = phi_I * ((I - I(-1)) - beta * (I(+1) - I));
    // NE8
    // q_n - pI = phi_I * ((I_n - I_n(-1)) - beta * (I_n(+1) - I_n));
    q_n = phi_I * ((I_n - I_n(-1)) - beta * (I_n(+1) - I_n));

    // LH5 
    K(+1) = (1 - delta) * K + delta * I;
    // NE7
    K_n(+1) = (1 - delta) * K_n + delta * I_n;

    /////////////////////////////////////////////////////////////////////////// 
    // --- FIRM BLOCK
    // Sticky equations start with LH, natural economy equations start with NE
    /////////////////////////////////////////////////////////////////////////// 
    // LF1 (redundant because we have an equation for the marginal cost in LF4)
    // Y = A + alpha * Omega + alpha * K + (1 - alpha) * N;

    // LF2
    rk = mc + Y - K;
    // NE5
    rk_n = Y_n - K_n;

    // LF3
    Wr = mc + Y - N;
    // NE4
    Wr_n = Y_n - N_n;

    // LF4
    // mc = -A + alpha * rk - alpha * Omega + (1 - alpha) * Wr;
    mc = alpha * rk - alpha * Omega + (1 - alpha) * Wr;
    // NE3
    // 0 = - A + alpha * rk_n - alpha * Omega + (1 - alpha) * Wr_n;
    0 = alpha * rk_n - alpha * Omega + (1 - alpha) * Wr_n;

    // LF5
    pi = beta * pi(+1) + (1 - theta) * (1 - theta * beta) / theta * mc;
    

    /////////////////////////////////////////////////////////////////////////// 
    // --- MONETARY POLICY
    /////////////////////////////////////////////////////////////////////////// 
    R = rho_R * R(-1) + (1 - rho_R) * (Phi_pi * pi + Phi_y * Y_g) + nu;
    nu = rho_nu * nu(-1) + eps_nu;

    /////////////////////////////////////////////////////////////////////////// 
    // --- EQUILIBRIUM 
    // Sticky equations start with LH, natural economy equations start with NE
    /////////////////////////////////////////////////////////////////////////// 
    // Y = s_c * c + s_I * I + s_I * pI;
    Y = s_c * c + s_I * I;
    // Y_n = s_c * c_n + s_I * I_n + s_I * pI;
    Y_n = s_c * c_n + s_I * I_n;

    /////////////////////////////////////////////////////////////////////////// 
    // --- GAP DEFINITIONS 
    /////////////////////////////////////////////////////////////////////////// 
    c_g = c - c_n;
    I_g = I - I_n;
    K_g = K - K_n; 
    Y_g = Y - Y_n;
    N_g = N - N_n; 
    Wr_g = Wr - Wr_n;
    q_g = q - q_n;
    rk_g = rk - rk_n;
end;

/*------------------------------------------------------------------------------
SHOCK VARIANCES / STANDARD DEVIATIONS
------------------------------------------------------------------------------*/
shocks;
    // With the TCA application already in mind, we will set the standard 
    // deviations of all except the monetary policy shock to a very small value
    var eps_z; stderr 1e-3;
    // var eps_pI; stderr 1e-3;
    // var eps_A; stderr 1e-3;
    var eps_Omega; stderr 1e-3;
    var eps_nu; stderr 1;
end;

/*------------------------------------------------------------------------------
MODEL CHECKS
------------------------------------------------------------------------------*/
resid;
steady;
check;

/*------------------------------------------------------------------------------
OBSERVED VARIABLES
Defined for TCA
------------------------------------------------------------------------------*/
varobs c I pi;

/*------------------------------------------------------------------------------
IRFs
------------------------------------------------------------------------------*/
// stoch_simul(order=1, irf=24) c_g I_g Y_g pi R nu Omega A pI z;
// stoch_simul(order=1, irf=24) c I Y pi R nu Omega A pI z;
// stoch_simul(order=1, irf=40, nograph);
stoch_simul(order=1, irf=12, nograph) R c I Keff pi nu Omega z;

