%% debugging get_aerosol_coefficients
clear all;

kappa = 0.6;
rd = 10.0;
T = 285;
P = 95000;

V = 1.0e-6; % 1 cm3

  cp = 1005               % Specific heat of air: J/kg K 
  Mw = 0.018              % Molecular weight of water: kg/mol
  Ma = 0.029              % Molecular weight of dry air: kg/mol
  g = 9.8                 % acceleration due to gravity: m/s^2
  R = 8.314               % Universal gas constant: J/K mol
  sigma_water = 0.072225  % surface tension of water (N/m)
  rho_w = 997             % density of water (kg/m3)

  P_atm = P/101325        % Pressure (atm)
  Dv = (0.211/P_atm) * (T/273)^(1.94)*(1e-4) % Mass diffusivity of water in air (m2/s or J/kg)
  
  % temperature-dependent parameters
  temp_c = T - 273.15
  a0 = 6.107799
  a1 = 4.436518e-1
  a2 = 1.428945e-2
  a3 = 2.650648e-4
  a4 = 3.031240e-6
  a5 = 2.034081e-8
  a6 = 6.136829e-11
  % vapor pressure of water (Pa)
  Po = 100*(a0+a1*temp_c+a2*(temp_c^2)+a3*(temp_c^3)+a4*(temp_c^4)+a5*(temp_c^5)+a6*(temp_c^6))
  % thermal conductivity of air (W/m K)                                          
  ka = 1e-3*(4.39+0.071*T)
  % density of air (kg/m3)
  rho_a = P/(287.058*T)
  % latent heat of vaporization: J/kg
  Hv = (2.5*((273.15/T)^(0.167+3.67e-4*T)))*1e6
  
  % Generalized coefficients
  G = 1/((rho_w*R*T/Po/Dv/Mw) + (Hv*rho_w/ka/T/(Hv*Mw/T/R - 1))) * 1e18     %nm2/sec
  A = 2*Mw*sigma_water/R/T/rho_w *1e9                                       %nm
  alpha = Hv*Mw*g/cp/R/T^2 - g*Ma/R/T                                       %1/m
  gamma1 = P*Ma/Po/Mw + Hv^2*Mw/cp/R/T^2
  gamma2 = gamma1*4*pi/rho_a/rho_w/V*1e-27                                   %1/nm3

  % 3-moment ODE coefficients
  a = zeros(5,1);
  a(1) = G;               %nm2/sec
  a(2) = G*A;             %nm3/sec
  a(3) = -G*kappa*rd^3;   %nm5/sec 
  a(4) = alpha;           %1/m
  a(5) = gamma2;          %1/nm3

  %%
  clear gamma
  N = 100.0;
  k = 3.8564;
  theta = 28.0704505;
  M = zeros(8,1);
  S = 0.000;
  v_up = 1;
  s = 5;
 
  for q=-4:2
      M(q+s) = N*theta^q*gamma(k+q)/gamma(k);
  end
  
  dMdt = zeros(3,1);
  
  for p=1:3
      dMdt(p+1) = p*(a(1)*M(p-2+s)*S+a(2)*M(p-3+s)+a(3)*M(p-5+s))
  end
  
  dSdt = alpha*v_up - gamma2*dMdt(end)/3