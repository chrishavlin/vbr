function VBR = CB_004_xfit_premelt()
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CB_004_xfit_premelt.m
  %
  %  Calls VBR using xfit_premelt method from:
  %    Hatsuki Yamauchi and Yasuko Takei, JGR 2016, "Polycrystal anelasticity at
  %    near-solidus temperatures,"
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %% write method list %%
  VBR.in.elastic.methods_list={'anharmonic','anh_poro'};
  VBR.in.anelastic.methods_list={'xfit_premelt'};

  % adjust Gu_0_ol to match YT2016
  VBR.in.elastic.anharmonic.Gu_0_ol=72.45; %[GPa]

  % create and use a custom anharmonic scaling to match derivatives
  % from YT2016. scaling for bulk modulus is arbitrary here, but
  % will not affect shear velocity and modulus calculations (unrelaxed
  % or relaxed).
  YT2016_derivatives.dG_dT = -10.94*1e6;
  YT2016_derivatives.dG_dP = 1.987; % Pa/Pa
  YT2016_derivatives.dG_dP2 = 0;

  YT2016_derivatives.dK_dT = 1.2 * YT2016_derivatives.dG_dT;
  YT2016_derivatives.dK_dP = 3 * YT2016_derivatives.dG_dP;
  YT2016_derivatives.dK_dP2 = 0;

  VBR.in.elastic.anharmonic.YT2016_derivatives = YT2016_derivatives;
  VBR.in.elastic.anharmonic.temperature_scaling = 'YT2016_derivatives';
  VBR.in.elastic.anharmonic.pressure_scaling = 'YT2016_derivatives';

  %% Define the Thermodynamic State %%

  VBR.in.SV.T_K=700:50:1200;
  VBR.in.SV.T_K=VBR.in.SV.T_K+273;
  sz=size(VBR.in.SV.T_K); % temperature [K]
  VBR.in.SV.P_GPa = 0.2 * ones(sz); % pressure [GPa]

  % this method requires the solidus.
  % you should write your own function for the solidus that takes all the other
  % state variables as input. This is just for illustration
  dTdz=0.5 ; % solidus slope [C/km]
  dTdP=dTdz / 3300 / 9.8 / 1000 * 1e9; % [C/GPa ]
  VBR.in.SV.Tsolidus_K=1000+dTdP*VBR.in.SV.P_GPa;

  % remaining state variables (ISV)
  VBR.in.SV.dg_um=3.1 * ones(sz); % grain size [um]
  VBR.in.SV.rho = 3300 * ones(sz); % density [kg m^-3]
  VBR.in.SV.sig_MPa = 10 * ones(sz); % differential stress [MPa]
  VBR.in.SV.phi = 0.0 * ones(sz); % melt fraction
  VBR.in.SV.f = 1./logspace(-2,4,100); % frequency range

  %% CALL THE VBR CALCULATOR %%
  [VBR] = VBR_spine(VBR) ;

  %% plot frequency dependence %%
  if ~vbr_tests_are_running()
    figure;
    subplot(1,3,1)
    semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.xfit_premelt.M(1,:,:)/1e9));
    ylabel('M [GPa]'); xlabel('period [s]')
    ylim([0,80])

    subplot(1,3,2)
    loglog(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.xfit_premelt.Qinv(1,:,:)));
    ylabel('Q^-1'); xlabel('period [s]')
    ylim([1e-3,.1])

    subplot(1,3,3)
    semilogx(1./VBR.in.SV.f,1e-3*squeeze(VBR.out.anelastic.xfit_premelt.V(1,:,:)));
    ylabel('V_s [km/s]'); xlabel('period [s]')
    saveas(gcf,'./figures/CB_004_xfit_premelt.png')
  end
end