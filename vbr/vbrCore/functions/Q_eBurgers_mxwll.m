function tau = Q_eBurgers_mxwll(VBR,Gu)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % tau = Q_eBurgers_mxwll(VBR,Gu)
  %
  % calculatues the maxwell time & limits for extended burgers model
  %
  % Parameters:
  % ----------
  % VBR    the VBR structure
  % Gu     unrelaxed modulus [GPa]
  %
  % Output:
  % ------
  % tau.   structure of maxwell times including:
  %    .maxwell = steady state viscous maxwell time (i.e., eta / Gu)
  %    .L = lower limit of integration for high temp background
  %    .H = upper limit of integration for high temp background
  %    .P = center period of dissipation peak (if being used)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % read in parameters
  Burger_params=VBR.in.anelastic.eburgers_psp;
  bType=Burger_params.eBurgerFit;
  fit_params = Burger_params.(bType);
  % state variables for either maxwell time or integration limits, peak loc:
  phi =  VBR.in.SV.phi ;
  T_K_mat = VBR.in.SV.T_K ; % temperature [K]
  P_Pa_mat = VBR.in.SV.P_GPa.*1e9 ; % convert pressure GPa to Pa = GPa*1e9
  d_mat = VBR.in.SV.dg_um ; % microns grain size
  Ch2o_ppm = VBR.in.SV.Ch2o;

  % Scaling values from JF10
  TR = fit_params.TR ;% Kelvins
  PR = fit_params.PR *1e9; % convert pressure GPa to Pa = GPa*1e9
  dR = fit_params.dR ; % microns grain size
  E = fit_params.E ; % activation energy J/mol
  if strcmp(bType, 'liu_water_2023')
      % the activation enthalpy depends on water conent
      E = E - fit_params.k * log(Ch2o_ppm);
  end
  R = Burger_params.R ; % gas constant
  Vstar = fit_params.Vstar ; % m^3/mol Activation Volume
  m_a = fit_params.m_a ; % grain size exponent (anelastic)
  m_v = fit_params.m_v ; % grain size exponent (viscous)

  % maxwell time calculation
  [visc_exists,missing]=checkStructForField(VBR,{'in','viscous','methods_list'},0);
  if Burger_params.useJF10visc || visc_exists==0
    % use JF10's exact relationship
    scale=((d_mat./dR).^m_v).*exp((E/R).*(1./T_K_mat-1/TR)).*exp((Vstar/R).*(P_Pa_mat./T_K_mat-PR/TR));
    scale=addMeltEffects(phi,scale,VBR.in.GlobalSettings,Burger_params);
    if strcmp(bType, 'liu_water_2023')
        scale = addWaterEffects(scale, Ch2o_ppm, fit_params);
    end
    Tau_MR = fit_params.Tau_MR ;
    tau.maxwell=Tau_MR .* scale ; % steady state viscous maxwell time
  else
    % use diffusion viscosity from VBR to get maxwell time
    visc_method=VBR.in.viscous.methods_list{1};
    eta_diff = VBR.out.viscous.(visc_method).diff.eta ; % viscosity for maxwell relaxation time
    tau.maxwell = eta_diff ./ Gu ; % maxwell relaxation time
  end

  % integration limits and peak location
  LHP=((d_mat./dR).^m_a).*exp((E/R).*(1./T_K_mat-1/TR)).*exp((Vstar/R).*(P_Pa_mat./T_K_mat-PR/TR));
  LHP=addMeltEffects(phi,LHP,VBR.in.GlobalSettings,Burger_params);

  % note on water effects: Liu et al find no dependence for the peak position
  % on water content, so here, we calculate tau.P as normal before adding the
  % water effects.
  tau.P = fit_params.Tau_PR * LHP;
  if strcmp(bType, 'liu_water_2023')
      LHP = addWaterEffects(LHP, Ch2o_ppm , fit_params);
  end
  tau.L = fit_params.Tau_LR * LHP;  
  tau.H = fit_params.Tau_HR * LHP;

end

function scaleMat = addWaterEffects(scaleMat, Ch2o_ppm, liu_params);
    c_ref = liu_params.c_ref;
    c_factor = (Ch2o_ppm / c_ref) .^ liu_params.r_m;
    scaleMat = scaleMat.*c_factor;
end

function scaleMat=addMeltEffects(phi,scaleMat,GlobalSettings,Burger_params)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % scaleMat=addMeltEffects(phi,scaleMat,GlobalSettings,Burger_params)
  %
  % adds on Melt Effects
  %
  % Parameters:
  % ----------
  % phi              melt fraction
  % scaleMat         the initial maxwell time matrix
  % GlobalSettings   global settings structure with melt_enhancement flag
  % Burger_params    the parameter structure for burgers model
  %
  % Output:
  % ------
  % scaleMat        the maxwell time matrix, adjusted for small melt effect
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % sharper response than creep at phi_c, but less sensitive at higher phi (for HTB only)
  alpha = Burger_params.melt_alpha ; % post-critical melt fraction dependence
  phi_c = Burger_params.phi_c ; % critical melt fraction
  x_phi_c = Burger_params.x_phi_c ;% melt enhancement factor

  % x_phi_c adjustment ("nominally melt free" to truly melt free)
  % x_phi_c will be 1 if melt_enhancement is turned off (set to 0)
  scaleMat = scaleMat.* x_phi_c ;

  % add melt effects
  [scale_mat_prime] = sr_melt_enhancement(phi,alpha,x_phi_c,phi_c) ;
  scaleMat = scaleMat ./ scale_mat_prime;
end
