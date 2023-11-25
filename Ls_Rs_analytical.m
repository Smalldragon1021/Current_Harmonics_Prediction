function [Ls,Rs] = Ls_Rs_analytical(Rs_0,Ls_0,leakage_factor,omega,d_Ls,ro_cu,mu0,airgap_factor_Ls,zt_Ls)


Lm_0 = (1-leakage_factor)*Ls_0;
Lsl_0 = leakage_factor*Ls_0;
% equivalent skin depth and reduced conductor height
alpha_Ls = sqrt(airgap_factor_Ls*omega*mu0/2/ro_cu); % inductance is proximity effect, resistance may be skin effect.
ksi_Ls = alpha_Ls.*d_Ls;
% inductance factors
phei_Ls = 3./(2.*ksi_Ls).*(sinh(2.*ksi_Ls)-sin(2.*ksi_Ls))./(cosh(2.*ksi_Ls)-cos(2.*ksi_Ls));
psi_Ls = 1./ksi_Ls.*(sinh(ksi_Ls)+sin(ksi_Ls))./(cosh(ksi_Ls)+cos(ksi_Ls));
k_L_Ls = 1/zt_Ls^2.*phei_Ls + (zt_Ls^2-1)/zt_Ls^2.*psi_Ls;
% resistance factors
fi_Ls = ksi_Ls.*(sinh(2.*ksi_Ls)+sin(2.*ksi_Ls))./(cosh(2.*ksi_Ls)-cos(2.*ksi_Ls));
sai_Ls = 2*ksi_Ls.*(sinh(ksi_Ls)-sin(ksi_Ls))./(cosh(ksi_Ls)+cos(ksi_Ls));
k_R_Ls = fi_Ls + (zt_Ls^2-1)/3.*sai_Ls;

Rs = Rs_0.*k_R_Ls;
Ls = Lm_0 + Lsl_0.*k_L_Ls;