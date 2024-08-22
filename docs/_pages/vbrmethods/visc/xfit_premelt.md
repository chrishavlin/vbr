---
permalink: /vbrmethods/visc/xfitpremelt/
title: ''
---

# xfit_premelt

Steady state flow law for pre-melting viscosity drop following Yamauchi and Takei (2016), J. Geophys. Res. Solid Earth, [DOI](https://doi.org/10.1002/2016JB013316). By default, this method uses flow law paramters (reference viscosity, activation energy, etc.) from a fit to the upper mantle.

## Requires

In addition to the state variable arrays required by the other viscous methods (see [HK2003](/vbr/vbrmethods/visc/hk2003/), [HZK2011](/vbr/vbrmethods/visc/hzk2011/)), this method also requires a solidus:

```matlab
VBR.in.SV.Tsolidus_K % solidus in degrees K
```

## Calling Procedure

```matlab
% set required state variables
clear
VBR.in.SV.T_K = linspace(800,1200,10)+273; % temperature [K]
VBR.in.SV.Tsolidus_K = (1200+273) * ones(size(VBR.in.SV.T_K)); % solidus [K]
VBR.in.SV.P_GPa = 2 * ones(size(VBR.in.SV.T_K)); % pressure [GPa]
VBR.in.SV.dg_um = 4 * ones(size(VBR.in.SV.T_K)); % grain size [um]
VBR.in.SV.sig_MPa = 10 *  ones(size(VBR.in.SV.T_K)); % differential stress [MPa]
VBR.in.SV.phi = 0.01 * ones(size(VBR.in.SV.T_K)); % melt fraction / porosity

% add to viscous methods list
VBR.in.viscous.methods_list={'xfit_premelt'};

% call VBR_spine
[VBR] = VBR_spine(VBR) ;
```

## Parameters

Parameters, including flow law parameters can be explored by loading them into the workspace with

```matlab
VBR.in.viscous.xfit_premelt = Params_Viscous('xfit_premelt');
disp(VBR.in.viscous.xfit_premelt)
```

### eta_melt_free_method

The parameter, `VBR.in.viscous.xfit_premelt.eta_melt_free_method`, controls what viscosity method is used to calculate 
the melt free viscosity. By default, this is set to `xfit_premelt`, in which case it uses the upper mantle fit directly 
from Yamauchi and Takei (2016). If set to one of the other viscosity mtehods, `HK2003` or `HZK2011`, then the melt 
free viscosity is calculated using those methods with melt fraction set to 0 and then the near-solidus pre-melting 
effect is then multiplied on.


### Grain size dependence 

Note that in section 4.4 of Yamauchi and Takei (2016), YT2016 fit for H, V following 
Priestly & McKenzie (2013) with `dg_um = dg_um_r` (i.e., the grain size is at the 
reference grain size). This results in a viscosity relationship independent of grain size. 
The VBRc sets `dg_um_r` to 4mm, the mean grain size of the upper mantle calculated in
Priestly & McKenzie (2013). To **exactly** match Yamauchi and Takei (2016), you should 
set the grain size to the value of the reference grain size (e.g., `VBR.in.SV.dg_um = 4000.`).

```matlab
VBR.in.viscous.xfit_premelt = Params_Viscous('xfit_premelt');
VBR.in.viscous.xfit_premelt.m = 0;
```

## Output
Output is stored in `VBR.out.viscous.xfit_premelt`. Unlike the other viscous methods, `xfit_premelt` only returns a diffusion creep viscosity sub-structure:

```matlab
>> disp(fieldnames(VBR.out.viscous.xfit_premelt))

{
  [1,1] = diff
}
>> disp(fieldnames(VBR.out.viscous.xfit_premelt.diff))

{
  [1,1] = eta
  [2,1] = eta_meltfree
}
```

Note that `eta_meltfree` is both the melt-free **and** volatile-free viscosity: the pre-melting method 
incorprates weakening from volatiles in the solidus depression. To incorporate effects of volatiles, 
should calculate a volatile-dependent solidus. 
