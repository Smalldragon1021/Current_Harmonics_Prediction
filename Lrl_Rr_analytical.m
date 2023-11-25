function [Lrl,Rr] = Lrl_Rr_analytical(Rr_0,Lrl_0,omega,d_Lrl,ro_cu,mu0,airgap_factor_Lrl)


% equivalent skin depth and reduced conductor height
alpha_Lrl = sqrt(airgap_factor_Lrl*omega*mu0/2/ro_cu); % inductance is proximity effect, resistance may be skin effect.
ksi_Lrl = alpha_Lrl.*d_Lrl;
% inductance factors
phei_Lrl = 3./(2.*ksi_Lrl).*(sinh(2.*ksi_Lrl)-sin(2.*ksi_Lrl))./(cosh(2.*ksi_Lrl)-cos(2.*ksi_Lrl));
%psi_Lrl = 1./ksi_Lrl.*(sinh(ksi_Lrl)+sin(ksi_Lrl))./(cosh(ksi_Lrl)+cos(ksi_Lrl));
k_L_Lrl = phei_Lrl;
% resistance factors
fi_Lrl = ksi_Lrl.*(sinh(2.*ksi_Lrl)+sin(2.*ksi_Lrl))./(cosh(2.*ksi_Lrl)-cos(2.*ksi_Lrl));
%sai_Lrl = 2*ksi_Lrl.*(sinh(ksi_Lrl)-sin(ksi_Lrl))./(cosh(ksi_Lrl)+cos(ksi_Lrl));
k_R_Lrl = fi_Lrl;

Rr = Rr_0.*k_R_Lrl;
Lrl =Lrl_0.*k_L_Lrl;