function [ VBR ] = wang2006( VBR )
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % [ VBR ] = wang2006( VBR )
  %
  % Conductivity of synthetic polycrystaline Ol (hydrous and Anhydrous)
  %
  % Parameters:
  % ----------
  % VBR    the VBR structure
  %
  % Output:
  % ------
  % VBR    the VBR structure, with VBR.out.electric.wang2006_ol structure
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % read in eletric parameters
  ele = VBR.in.electric.wang2006_ol;
  T = VBR.in.SV.T; % K (Temperature)
  Ch2o = VBR.in.SV.Ch2o; % ppm (water content)
  P = VBR.in.SV.P; % GPa (Pressure)
  
    % Hydrous Conduction
    SH = ele.S_H; % S/m
    Va_H = ele.Va_H; % cc/mol
    HH = ele.H_H + Va_H.*P; % kJ
    R = ele.R_H; % kJ/(mol*K)
    a = ele.a_H; % unitless
    r = ele.r_H; % unitless
    
    % Anhydrous Conduction
    SA = ele.S_A; % S/m
    Va_A = ele.Va_A; % cc/mol
    HA = ele.H_A + Va_A.*P; % kJ

  % calculate arrhenius relation for each conduction mechanism
  esig_A = arrh_dry(SA,HA,R,T);
  esig_H = arrh_wet(SH,HH,R,T,Ch2o,a,r);
  
  % summation of conduction mechanisms
  esig = esig_A + esig_H; % S/m
  
  % store in VBR structure
  wang2006_ol.esig_i = esig_i;
  wang2006_ol.esig_h = esig_h;
  wang2006_ol.esig_p = esig_p;
  wang2006_ol.esig = esig;
  VBR.out.electric.wang2006_ol = wang2006_ol;
end

function sig = arrh_dry(S,H,k,T)
    exponent = -(H)./(k.*T);
    sig = (10^S).*exp(exponent);
end

function sig = arrh_wet(S,H,k,T,w,a,r)
 exponent = -(H-a.*(w.^(1/3)))./(k.*T);
    sig = (10^S).*(w.^r).*exp(exponent);
end