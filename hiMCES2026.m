%===============================================================================
% Version JMCC PLUS (DOI: XXX)
% Author: Amin Forouzandehmehr
%===============================================================================

function [dY, dati] = hiMCES2026(time, Y, contrflag, GAC, stimFlag, freq, classicVSoptic, tDrugApplication, INaFRedMed, INaLRedMed, INaKRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, IfRedMed, SERCARedMed, SERCAKmRedMed, RyRRedMed,RyRPoRedMed, ItoRedMed, ischemiaFlag, myFlags, IKatpSelectionString, fnp, z, flv, fmyo, gma, ISOFlag)

%-------------------------------------------------------------------------------
% State variables
%-------------------------------------------------------------------------------

% Electrophysiology:

% 1: Vm (volt) (in Membrane)
% 2: Ca_SR (millimolar) (in calcium_dynamics)
% 3: Cai (millimolar) (in calcium_dynamics)
% NOT USED 4: g (dimensionless) (in calcium_dynamics)
% 5: d (dimensionless) (in i_CaL_d_gate)
% 6: f1 (dimensionless) (in i_CaL_f1_gate)
% 7: f2 (dimensionless) (in i_CaL_f2_gate)
% 8: fCa (dimensionless) (in i_CaL_fCa_gate)
% 9: Xr1 (dimensionless) (in i_Kr_Xr1_gate)
% 10: Xr2 (dimensionless) (in i_Kr_Xr2_gate)
% 11: Xs (dimensionless) (in i_Ks_Xs_gate)
% 12: h (dimensionless) (in i_Na_h_gate)
% 13: j (dimensionless) (in i_Na_j_gate)
% 14: m (dimensionless) (in i_Na_m_gate)
% 15: Xf (dimensionless) (in i_f_Xf_gate)
% 16: q (dimensionless) (in i_to_q_gate)
% 17: r (dimensionless) (in i_to_r_gate)
% 18: Nai (millimolar) (in sodium_dynamics)
% 19: m_L (dimensionless) (in i_NaL_m_gate)
% 20: h_L (dimensionless) (in i_NaL_h_gate)
% 21: RyRa (dimensionless) (in calcium_dynamics)
% 22: RyRo (dimensionless) (in calcium_dynamics)
% 23: RyRc (dimensionless) (in calcium_dynamics)

% Myofilament:

% 24: Nxb         Non-permissive fraction in overlap region
% 25: XBpreR      Fraction of XBs in pre-rotated state (Stiffness generation)
% 26: XBpostR     Fraction of XBS in post-rotated state (force generation)
% 27: x_XBpreR    Strain of XBs in the pre-rotated state
% 28: x_XBpostR   Strain of XBs in the post-rotated state
% 29: TropCaL     Ca bound to the low affinity troponin site
% 30: TropCaH     Ca bound to the high affinity troponin site
% 31: IntegF      Integral of force
% 32: SL          Sarcomere length

% Oxygen dynamics
% 33: O2o         Extracellular O2 concentration

%% Constants
F = 96485.3415;     % coulomb_per_mole (in model_parameters)
R = 8.314472;       % joule_per_mole_kelvin (in model_parameters)
if classicVSoptic == 0
    T = 310.0;      % kelvin (in model_parameters) %37°C
elseif classicVSoptic == 1
    T = 294.0;      % kelvin (in model_parameters) %21°C
end

%% Q10 factors 
expDeltaT = (310-T)/10;
Q10_INa         = 2;       coeff_Q10_INa_tau       = Q10_INa^expDeltaT;
Q10_INaL        = 2.2;     coeff_Q10_INaL_tau      = Q10_INaL^expDeltaT;
Q10_ICaL        = 2.1;     coeff_Q10_ICaL_tau      = Q10_ICaL^expDeltaT;
Q10_If          = 4.5;     coeff_Q10_If_tau        = Q10_If^expDeltaT;
Q10_IKr_act     = 4.55;    coeff_Q10_IKr_act_tau   = Q10_IKr_act^expDeltaT;     %Mauerhofer x Q10 IKr
Q10_IKr_inact   = 3.08;    coeff_Q10_IKr_inact_tau = Q10_IKr_inact^expDeltaT;   %media 2.94 e 3.22
Q10_generico    = 2;       coeff_Q10_generico_tau  = Q10_generico^expDeltaT;

%% Cell geometry
V_SR = 583.73;        % micrometre_cube (in model_parameters)
Vc   = 8800.0;        % micrometre_cube (in model_parameters)
Cm   = 9.87109e-11;   % farad (in model_parameters)

%% Extracellular concentrations
if classicVSoptic==0
    % ORIGINAL
    Nao = 151.0; % millimolar (in model_parameters)
    Ko  = 5.4;   % millimolar (in model_parameters)    
    Cao = 1.8;   % millimolar (in model_parameters)
elseif classicVSoptic == 1
    % Optical measurements
    Nao = 135.0;  % millimolar (in model_parameters)
    Ko  = 5.4;    % millimolar (in model_parameters)
    Cao = 1.33;   % millimolar (in model_parameters)
end

%% Ischemia
% All parameters and references: Table 1 in the manuscript.

flev = flv; % fKATP for different concentrations of Levosimendan. 10 uM Levo = 1.843903, 2 uM = 1.118538, 0.3 uM = 1.002385. 
flev = 1;   % No Levosimendan simulation
if time >= 350 && time < 850

if ischemiaFlag == 0
    % Normoxia
    fkatp = 0.0041 * ((time<tDrugApplication)*1+(time >= tDrugApplication)*flev);
    finhib = 1;
    fnk = fnp; fpc = fnp;
    myFlags = zeros(1,7);   % regardless of user's inputs

elseif ischemiaFlag == 1    %***********************************************************************************************************
    % Mild 
    
    %Hyperkalemia OFF:0 / ON:1
    if myFlags(1)==1
        Ko =6.25;      %mM         Ko = 6.25
    end
    
    %Acidosis OFF:0 / ON:1
    if myFlags(2)==1
        finhib = 0.875; %Ledezma - Zhou 10.1371/journal.pone.0220294
    elseif myFlags(2)==0
        finhib = 1;
    end
    
    %Hypoxia IKatp OFF:0 / ON:1
    if myFlags(3)==1
        if isequal(IKatpSelectionString, 'Ledezma')
            fkatp = 0.0092; %0.1;
        elseif isequal(IKatpSelectionString, 'Kazbanov')
            fkatp = 0.0092 * ((time<tDrugApplication)*1+(time >= tDrugApplication)*flev); % 0.0055/2  * ((time<tDrugApplication)*1+(time >= tDrugApplication)*5); % half of the severity 2 fkatp by Kazbanov, they didn't have an intermediate value. Same ratio as in Ledezma. *6 for Levo
        end
    elseif myFlags(3)==0
        fkatp = 0.0041;
    end
    
    
    %Hypoxia on INaK OFF:0 / ON:1
    if myFlags(4)==1
        fnk = fnp;

    elseif myFlags(4)==0
          fnk = 0;
    end
    
    %Hypoxia on IpCa OFF:0 / ON:1
    if myFlags(5)==1
        fpc = fnp;

    elseif myFlags(5)==0
          fpc = 0;
    end
    
elseif ischemiaFlag == 2    % **************************************************************************************************************
    % Severe
    
    %Hyperkalemia OFF:0 / ON:1
    if myFlags(1)==1
    Ko = 9;         %mM
    end
        
    %Acidosis OFF:0 / ON:1
    if myFlags(2)==1
        finhib = 0.75;  %Ledezma - Zhou 10.1371/journal.pone.0220294
    elseif myFlags(2)==0
        finhib = 1;
    end
    
    %Hypoxia IKatp OFF:0 / ON:1
    if myFlags(3)==1
        if isequal(IKatpSelectionString, 'Ledezma')
            fkatp = 0.0128; %0.2;
        elseif isequal(IKatpSelectionString, 'Kazbanov')
            fkatp = 0.0128 * ((time<tDrugApplication)*1+(time >= tDrugApplication)*flev); %0.0055; Levosimendan effect on IfKATP
        end
    elseif myFlags(3)==0
        fkatp = 0.0041;
    end

    %Hypoxia on INaK OFF:0 / ON:1
    if myFlags(4)==1
        fnk = fnp;

    elseif myFlags(4)==0
          fnk = 0;
    end
    
    %Hypoxia on IpCa OFF:0 / ON:1
    if myFlags(5)==1
        fpc = fnp;

    elseif myFlags(5)==0
          fpc = 0;
    end
    
end
else if time >= 850     % Reoxygenation
    % Normoxia
    fkatp = 0.0041 * ((time<tDrugApplication)*1+(time >= tDrugApplication)*flev);
    finhib = 1;
    fnk = fnp; fpc = fnp;
    myFlags = zeros(1,7);   
end
end

if time < 350
        % Normoxia
    fkatp = 0.0041 * ((time<tDrugApplication)*1+(time >= tDrugApplication)*flev);
    finhib = 1;
    fnk = fnp; fpc = fnp;
    myFlags = zeros(1,7);   %disregarding user's input
end

%% Intracellular concentrations
% Naio = 10 mM Y(18)
Ki = 150.0;   % millimolar (in model_parameters)
% Cai  = 0.0002 mM Y(3)
% caSR = 0.3 mM Y(2)
% time (second)

%% Nernst potential
E_Na = R*T/F*log(Nao/Y(18));
E_Ca = 0.5*R*T/F*log(Cao/Y(3));
E_K  = R*T/F*log(Ko/Ki);
PkNa = 0.03;   % dimensionless (in electric_potentials)
E_Ks = R*T/F*log((Ko+PkNa*Nao)/(Ki+PkNa*Y(18)));

%% INa 
g_Na        = GAC(1)* ((time<tDrugApplication)*1+(time >= tDrugApplication)*INaFRedMed)*finhib*6447.1896;
i_Na        =  g_Na*Y(14)^3.0*Y(12)*Y(13)*(Y(1) - E_Na);

m_inf       = 1 / (1 + exp((Y(1)*1000 + 39)/-11.2));
tau_m       = GAC(17) * coeff_Q10_INa_tau* (0.00001 + 0.00013*exp(-((Y(1)*1000 + 48)/15)^2) + 0.000045 / (1 + exp((Y(1)*1000 + 42)/-5)));
dY(14, 1)   = (m_inf-Y(14))/tau_m;

h_inf       = 1 / (1 + exp((Y(1)*1000 + 66.5)/6.8));
tau_h       = GAC(18) * coeff_Q10_INa_tau * (0.00007 + 0.034 / (1 + exp((Y(1)*1000 + 41)/5.5) + exp(-(Y(1)*1000 + 41)/14)) + 0.0002 / (1 + exp(-(Y(1)*1000 + 79)/14)));
dY(12, 1)   = (h_inf-Y(12))/tau_h;

j_inf       = h_inf;
tau_j       = GAC(18) * coeff_Q10_INa_tau * 10*(0.0007 + 0.15 / (1 + exp((Y(1)*1000 + 41)/5.5) + exp(-(Y(1)*1000 + 41)/14)) + 0.002 / (1 + exp(-(Y(1)*1000 + 79)/14)));
dY(13, 1)   = (j_inf-Y(13))/tau_j;


%% INaL
myCoefTauM  = 1;
tauINaL     = 200; %ms
GNaLmax     = GAC(11) * ((time<tDrugApplication)*1+(time >= tDrugApplication)*INaLRedMed)* finhib*2.3*7.5; %(S/F)
Vh_hLate    = 87.61;
i_NaL       = GNaLmax* Y(19)^(3)*Y(20)*(Y(1)-E_Na);

m_inf_L     = 1/(1+exp(-(Y(1)*1000+42.85)/(5.264)));
alpha_m_L   = 1/(1+exp((-60-Y(1)*1000)/5));
beta_m_L    = 0.1/(1+exp((Y(1)*1000+35)/5))+0.1/(1+exp((Y(1)*1000-50)/200));
tau_m_L     = coeff_Q10_INaL_tau * 1/1000 * myCoefTauM*alpha_m_L*beta_m_L;
dY(19, 1)   = (m_inf_L-Y(19))/tau_m_L;

h_inf_L     = 1/(1+exp((Y(1)*1000+Vh_hLate)/(7.488)));
tau_h_L     = coeff_Q10_INaL_tau * 1/1000 * tauINaL;
dY(20, 1)   = (h_inf_L-Y(20))/tau_h_L;

%% If 

ISO_on        = time >= tDrugApplication && ISOFlag == 1;   % step; replace with a smooth ramp if desired
dVhalf_mV     = 8;      % +5 to +12 mV typical positive shift under ISO      ISO: 8
tau_scale_ISO = 1.0;    % <1 speeds gating                                   ISO: 1.0
ISO_level = double(ISO_on);

g_f         = GAC(2) * ((time<tDrugApplication)*1+(time >= tDrugApplication)*IfRedMed)*22.2763088;

fNa         = 0.37;
fK          = 1 - fNa;
i_fK        = fK*g_f*Y(15)*(Y(1) - E_K);
i_fNa       = fNa*g_f*Y(15)*(Y(1) - E_Na);

i_f         = i_fK + i_fNa;

Vhalf_base_mV = -69;                 
k_mV          = 8;                   % slope factor
V_mV          = Y(1)*1000;

Vhalf_mV = Vhalf_base_mV + dVhalf_mV * ISO_level;     % depolarizing shift
Xf_infinity = 1.0 ./ (1.0 + exp((V_mV - Vhalf_mV)/k_mV));

%Xf_infinity = 1.0/(1.0 + exp((Y(1)*1000 + 69)/8));
%tau_Xf      = coeff_Q10_If_tau * 5600 / (1 + exp((Y(1)*1000 + 65)/7) + exp(-(Y(1)*1000 + 65)/19));
tau_base = coeff_Q10_If_tau * 5600 ./ (1 + exp((V_mV + 65)/7) + exp(-(V_mV + 65)/19));
tau_Xf   = tau_base .* (1 - ISO_level + ISO_level * tau_scale_ISO);   % <= reduces τ under ISO


dY(15, 1)   = 1000*(Xf_infinity-Y(15))/tau_Xf;

%% ICaL
g_CaL       = GAC(3)* finhib*8.635702e-5;   % metre_cube_per_F_per_s (in i_CaL)
i_CaL       = ((time<tDrugApplication)*1+(time >= tDrugApplication)*ICaLRedMed)*g_CaL*4.0*Y(1)*F^2.0/(R*T)*(Y(3)*exp(2.0*Y(1)*F/(R*T))-0.341*Cao)/(exp(2.0*Y(1)*F/(R*T))-1.0)*Y(5)*Y(6)*Y(7)*Y(8);

d_infinity  = 1.0/(1.0+exp(-(Y(1)*1000.0+9.1)/7.0));
alpha_d     = 0.25+1.4/(1.0+exp((-Y(1)*1000.0-35.0)/13.0));
beta_d      = 1.4/(1.0+exp((Y(1)*1000.0+5.0)/5.0));
gamma_d     = 1.0/(1.0+exp((-Y(1)*1000.0+50.0)/20.0));
tau_d       = GAC(19) * coeff_Q10_ICaL_tau * (alpha_d*beta_d+gamma_d)*1.0/1000.0;
dY(5, 1)    = (d_infinity-Y(5))/tau_d;

f1_inf      = 1.0/(1.0+exp((Y(1)*1000.0+26.0)/3.0));
if (f1_inf-Y(6) > 0.0)
    constf1 = 1.0+1433.0*(Y(3)-50.0*1.0e-6);
else
    constf1 = 1.0;
end
tau_f1      = GAC(20) * coeff_Q10_ICaL_tau * (20.0+1102.5*exp(-((Y(1)*1000.0+27.0)/15.0)^2.0)+200.0/(1.0+exp((13.0-Y(1)*1000.0)/10.0))+180.0/(1.0+exp((30.0+Y(1)*1000.0)/10.0)))*constf1/1000.0;
dY(6, 1)    = (f1_inf-Y(6))/tau_f1;

f2_inf      = 0.33+0.67/(1.0+exp((Y(1)*1000.0+32.0)/4.0));
constf2     = 1.0;
tau_f2      = GAC(20) * coeff_Q10_ICaL_tau * (600.0*exp(-(Y(1)*1000.0+25.0)^2.0/170.0)+31.0/(1.0+exp((25.0-Y(1)*1000.0)/10.0))+16.0/(1.0+exp((30.0+Y(1)*1000.0)/10.0)))*constf2/1000.0;
dY(7, 1)    = (f2_inf-Y(7))/tau_f2;

alpha_fCa   = 1.0/(1.0+(Y(3)/0.0006)^8.0);
beta_fCa    = 0.1/(1.0+exp((Y(3)-0.0009)/0.0001));
gamma_fCa   = 0.3/(1.0+exp((Y(3)-0.00075)/0.0008));
fCa_inf     = (alpha_fCa+beta_fCa+gamma_fCa)/1.3156;
if ((Y(1) > -0.06) && (fCa_inf > Y(8)))
    constfCa = 0.0;
else
    constfCa = 1.0;
end
tau_fCa     = coeff_Q10_ICaL_tau * 0.002;   % second (in i_CaL_fCa_gate)
dY(8, 1)    = constfCa*(fCa_inf-Y(8))/tau_fCa;

%% Ito
g_to        = GAC(4) * ((time<tDrugApplication)*1+(time >= tDrugApplication)*ItoRedMed) * 29.9038;   % S_per_F (in i_to)
i_to        = g_to*(Y(1)-E_K)*Y(16)*Y(17);

q_inf       = 1.0/(1.0+exp((Y(1)*1000.0+53.0)/13.0));
tau_q       = coeff_Q10_generico_tau * (6.06+39.102/(0.57*exp(-0.08*(Y(1)*1000.0+44.0))+0.065*exp(0.1*(Y(1)*1000.0+45.93))))/1000.0;
dY(16, 1)   = (q_inf-Y(16))/tau_q;

r_inf       = 1.0/(1.0+exp(-(Y(1)*1000.0-22.3)/18.75));
tau_r       = coeff_Q10_generico_tau * (2.75352+14.40516/(1.037*exp(0.09*(Y(1)*1000.0+30.61))+0.369*exp(-0.12*(Y(1)*1000.0+23.84))))/1000.0;
dY(17, 1)   = (r_inf-Y(17))/tau_r;

%% IKs
%if time >= 350 && time < 850 && any (GAC(1:12)-1)
%ro2 =  0.15;
%else
    ro2 = 1;
%end

g_Ks        = GAC(5) * 2.041;   % S_per_F (in i_Ks)
i_Ks        = ro2* ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKsRedMed)*g_Ks*(Y(1)-E_Ks)*Y(11)^2.0*(1.0+0.6/(1.0+(3.8*0.00001/Y(3))^1.4));

Xs_infinity = 1.0/(1.0+exp((-Y(1)*1000.0-20.0)/16.0));
alpha_Xs    = 1100.0/sqrt(1.0+exp((-10.0-Y(1)*1000.0)/6.0));
beta_Xs     = 1.0/(1.0+exp((-60.0+Y(1)*1000.0)/20.0));
tau_Xs      = coeff_Q10_generico_tau * 1.0*alpha_Xs*beta_Xs/1000.0;
dY(11, 1)   = (Xs_infinity-Y(11))/tau_Xs;

%% IKr

%if time >= 350 && time < 850 && any (GAC(1:12)-1)
%ro1 =  0.15;
%else
    ro1 = 1;
%end

L0           = 0.025;   % dimensionless (in i_Kr_Xr1_gate)
Q            = 2.3;     % dimensionless (in i_Kr_Xr1_gate)
g_Kr         = GAC(6) * 29.8667;   % S_per_F (in i_Kr)
i_Kr         = ro1* ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKrRedMed)*g_Kr*(Y(1)-E_K)*Y(9)*Y(10)*sqrt(Ko/5.4);

V_half       = 1000.0*(-R*T/(F*Q)*log((1.0+Cao/2.6)^4.0/(L0*(1.0+Cao/0.58)^4.0))-0.019);

Xr1_inf      = 1.0/(1.0+exp((V_half-Y(1)*1000.0)/4.9));
alpha_Xr1    = 450.0/(1.0+exp((-45.0-Y(1)*1000.0)/10.0));
beta_Xr1     = 6.0/(1.0+exp((30.0+Y(1)*1000.0)/11.5));
tau_Xr1      = coeff_Q10_IKr_act_tau * 1.0*alpha_Xr1*beta_Xr1/1000.0;
dY(9, 1)     = (Xr1_inf-Y(9))/tau_Xr1;

Xr2_infinity = 1.0/(1.0+exp((Y(1)*1000.0+88.0)/50.0));
alpha_Xr2    = 3.0/(1.0+exp((-60.0-Y(1)*1000.0)/20.0));
beta_Xr2     = 1.12/(1.0+exp((-60.0+Y(1)*1000.0)/20.0));
tau_Xr2      = coeff_Q10_IKr_inact_tau * 1.0*alpha_Xr2*beta_Xr2/1000.0;
dY(10, 1)    = (Xr2_infinity-Y(10))/tau_Xr2;

%% IK1

alpha_K1    = 3.91/(1.0+exp(0.5942*(Y(1)*1000.0-E_K*1000.0-200.0)));
beta_K1     = (-1.509*exp(0.0002*(Y(1)*1000.0-E_K*1000.0+100.0))+exp(0.5886*(Y(1)*1000.0-E_K*1000.0-10.0)))/(1.0+exp(0.4547*(Y(1)*1000.0-E_K*1000.0)));
XK1_inf     = alpha_K1/(alpha_K1+beta_K1);
g_K1        = GAC(7) * 28.1492;   % S_per_F (in i_K1)
i_K1        = 1* g_K1*XK1_inf*(Y(1)-E_K)*sqrt(Ko/5.4);

%% INCX

KmCa        = 1.38;   % millimolar (in i_NaCa)
KmNai       = 87.5;   % millimolar (in i_NaCa)
Ksat        = 0.1;    % dimensionless (in i_NaCa)
gamma       = 0.35;   % dimensionless (in i_NaCa)
alpha       = 2.16659;
kNaCa      = ((time<tDrugApplication)*1+(time >= tDrugApplication)*INaCaRedMed) * GAC(8) * 6514.47574;   % A_per_F (in i_NaCa)
i_NaCa      = kNaCa*(exp(gamma*Y(1)*F/(R*T))*Y(18)^3.0*Cao-exp((gamma-1.0)*Y(1)*F/(R*T))*Nao^3.0*Y(3)*alpha)/((KmNai^3.0+Nao^3.0)*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*Y(1)*F/(R*T))));

%% INaK

ronak = 1-fnk*(0.5*(1+tanh((0.06)*(time-300))))+fnk*(0.5*(1+tanh((0.06)*(time-900)))); % 0.2 for SEV1 & 0.31 for SEV2

Km_K        = 1.0;    % millimolar (in i_NaK)
Km_Na       = ((time<tDrugApplication)*1+(time >= tDrugApplication)*INaKRedMed) * 40.0;   % millimolar (in i_NaK)
PNaK        = 2.74240 * GAC(9); % A_per_F (in i_NaK)
i_NaK       = 1 * ronak * PNaK*Ko/(Ko+Km_K)*Y(18)/(Y(18)+Km_Na)/(1.0+0.1245*exp(-0.1*Y(1)*F/(R*T))+0.0353*exp(-Y(1)*F/(R*T)));

%% IpCa

ropca = 1-fpc*(0.5*(1+tanh((0.06)*(time-300))))+fpc*(0.5*(1+tanh((0.06)*(time-900)))); % 0.2 for SEV1 & 0.31 for SEV2

KPCa        = 0.0005;   % millimolar (in i_PCa)

g_PCa       = GAC(10) * 0.4125;   % A_per_F (in i_PCa)
i_PCa       = ropca*g_PCa*Y(3)/(Y(3)+KPCa);

%% Background currents
g_b_Na      = 1.14;         % S_per_F (in i_b_Na)
i_b_Na      = g_b_Na*(Y(1)-E_Na);

g_b_Ca      = 0.8727264;    % S_per_F (in i_b_Ca)
i_b_Ca      = g_b_Ca*(Y(1)-E_Ca);

%% IKATP
% An ATP-sensitive K+ current (IKATP) was introduced to the hiMCES model based on a formulation by Kazbanov et al. (2014). The maximum IKATP conductance was rescaled by fKATP values in the control, SEV1 and SEV2 conditions (Table 1 of the paper). 
% Values of fKATP were calculated based on an ATP- and ADP-dependent approach by Ferrero et al. (1996).
K_o_n = 5.4; 
v = Y(1);
EK = E_K;

%fkatp = 0.0041;      % fkatp = 0.0041 for control simulations

if isequal(IKatpSelectionString, 'Ledezma')
    gkatp = 0.064;   % milliS_per_microF (in I_katp)
    akik = (Ko/K_o_n)^0.24;
    IKatp = fkatp*gkatp*akik*(v-EK);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IKatp Kazbanov 10.1371/journal.pcbi.1003891
if isequal(IKatpSelectionString, 'Kazbanov')
    gkatp = 155;   % milliS_per_microF (in I_katp)
    akik = (Ko/K_o_n)^0.3;
    bkik = 1/(40+3.5*exp(0.025*v));
    IKatp = fkatp*gkatp*akik*bkik*(v-EK);
end

%% Sarcoplasmic reticulum

% Metabolite-sensitive SERCA pump adapted from Tran et al., 2009
% (https://doi.org/10.1016/j.bpj.2008.11.045).

n = 1;
ATP = 5;
ADP = ones(1,n)*36e-3;
Pi = ones(1,n)*2;
Ca_i = Y(3);
Ca_sr = Y(2);
H_i = ones(1,n)*10^(-7.15)*1e3;    % Def: 7.15
sca = ones(1,8);
sck = ones(1,6);

    if ischemiaFlag == 1  && myFlags(6)==1 % SEV1 (5 min) % Hypoxia on Jup OFF:0 / ON:1
        ATP = 4.2;
        ADP = ones(1,n)*0.33;
        Pi = ones(1,n) * 15.5;
        H_i = ones(1,n)*10^(-6.5)*1e3;
        sca = [0.1, 1, 1, 1, 0.1, 5, 0.97, 1]; %kdcai, kdcasr, kdh1, kdhi, kdhsr, kdh, n_h, kdatp
        sck = [2  2  5  1  0.1  1];            %sck = [2  2  5  1  0.1  1]; kp1, kp2, kp3, km1=1, km2, km3
        % For param6 sca(1) = 0.1
        % For param6 sca(7) = 0.88
        
    elseif ischemiaFlag == 2 && myFlags(6)==1   % SEV2 (10 min) % Hypoxia on Jup OFF:0 / ON:1
        ATP = 3.1;
        ADP = ones(1,n)*0.25;
        Pi = ones(1,n)*22;
        H_i = ones(1,n)*10^(-6.3)*1e3;
        sca = [0.1, 1, 1, 1, 0.1, 5, 0.97, 1]; % kdcai, kdcasr, kdh1, kdhi, kdhsr, kdh, n_h, kdatp,
        sck = [2  2  5  1  0.1  1];      %sck = [2  2  5  1  0.1  1]; kp1, kp2, kp3, km1=1, km2, km3
    end

        % Size of input vectors
        n = size(ATP,2);

        %Creating a vector of the parameters
        k_p1 = 25900*ones(1,n) * sck(1);
        k_p2 = 2540*ones(1,n) * sck(2);
        k_p3 = 20.5*ones(1,n) * sck(3);
        k_m1 = 2*ones(1,n) * sck(4);
        k_m2 = 67200*ones(1,n) * sck(5);
        k_m3 = 149*ones(1,n) * sck(6);
        
        kdcai = 0.91*ones(1,n) * sca(1) * ((time<tDrugApplication)*1+(time >= tDrugApplication)*SERCAKmRedMed) * GAC(22);
        kdcasr = 2.24*ones(1,n) * sca(2);
        kdh1 = 1.09e-5*ones(1,n) * sca(3);
        kdhi = 3.54e-3*ones(1,n) * sca(4);
        kdhsr = 1.05e-8*ones(1,n) * sca(5);
        kdh = 7.24e-5*ones(1,n)  * sca(6);
        n_H = 2 * sca(7);   % Def: 2
        
        kdatp = (k_m1./k_p1) * sca(8);
        
        T_ATP = ATP./kdatp;
        T_Cai = Ca_i./kdcai;
        T_Casr = Ca_sr./kdcasr;
        T_H1 = H_i./kdh1;
        T_Hi = (H_i.^n_H)./kdhi;
        T_Hsr = (H_i.^n_H)./kdhsr;
        T_H = H_i./kdh;

        % Denominator terms for the apparent rate constants
        denom1 = (T_Hi + T_ATP.*T_Hi + T_ATP.*T_Hi.*T_H1 + T_ATP.*T_Cai.^2.*T_Hi + T_ATP.*T_Cai.^2);
        denom2 = (T_Casr.*T_Casr.*T_H + T_H.*(1) + T_Hsr.*(1+T_H));
        
        % Apparent forward rate constants 
        a_p1 = (k_p2.*T_ATP.*T_Cai.^2)./denom1;
        a_p2 = (k_p3.*T_Hsr)./denom2;
        
        % Apparent backward rate constants
        a_m1 = (k_m2.*ADP.*T_Casr.^2.*T_H)./denom2;
        a_m2 = (k_m3.*Pi.*T_Hi)./denom1;

        
        % Calculating steady state flux using the diagram method
        s1 = a_m1 + a_p2;
        s2 = a_p1 + a_m2;

        sum = s1 + s2;

        % Calculating the cycle rate
        v_cycle = (a_p1.*a_p2 - a_m1.*a_m2)./sum;

        % Proportion of states that are phosphorylated
        phosph = s2./sum;


i_up        = ((time<tDrugApplication)*1+(time >= tDrugApplication)*SERCARedMed) * GAC(15) * v_cycle * 0.0308;      % Def: 2*v_cycle*0.0154, for normal single cell: 0.0308.


V_leak		= 4.48209e-4;
i_leak      = (Y(2)-Y(3))*V_leak;

dY(4, 1)    = 0;

% RyR
g_irel_max	= 55.808061 * GAC(16) * ((time<tDrugApplication)*1+(time >= tDrugApplication)*RyRRedMed);
RyRa1       = 0.05169;
RyRa2       = 0.050001;
RyRahalf    = 0.02632;
RyRohalf    = 0.00944 * ((time<tDrugApplication)*1+(time >= tDrugApplication)*RyRPoRedMed);
RyRchalf    = 0.00167* GAC(21);

RyRSRCass   = (1 - 1/(1 +  exp((Y(2)-0.3)/0.1)));
i_rel       = g_irel_max*RyRSRCass*Y(22)*Y(23)*(Y(2)-Y(3));

RyRainfss   = RyRa1-RyRa2/(1 + exp((1000*Y(3)-(RyRahalf))/0.0082));
RyRtauadapt = GAC(12) * coeff_Q10_generico_tau * 1; %s
dY(21,1)    = (RyRainfss- Y(21))/RyRtauadapt;

RyRoinfss   = (1 - 1/(1 +  exp((1000*Y(3)-(Y(21)+ RyRohalf))/0.003)));
if (RyRoinfss>= Y(22))
    RyRtauact = GAC(13) * coeff_Q10_generico_tau * 18.75e-3;       %s
else
    RyRtauact = GAC(13) * coeff_Q10_generico_tau * 0.1*18.75e-3;   %s
end
dY(22,1)    = (RyRoinfss- Y(22))/(RyRtauact);

RyRcinfss   = (1/(1 + exp((1000*Y(3)-(Y(21)+RyRchalf))/0.001)));
if (RyRcinfss>= Y(23))
    RyRtauinact = GAC(14) * coeff_Q10_generico_tau * 2*87.5e-3;    %s
else
    RyRtauinact = GAC(14) * coeff_Q10_generico_tau * 87.5e-3;      %s
end
dY(23,1)    = (RyRcinfss- Y(23))/(RyRtauinact);

%% Ca2+ buffering
Buf_C       = 0.25-(70e-3);   % millimolar (in calcium_dynamics)
Buf_SR      = 10.0;   % millimolar (in calcium_dynamics)
Kbuf_C      = 0.001;   % millimolar (in calcium_dynamics)
Kbuf_SR     = 0.3;   % millimolar (in calcium_dynamics)
Cai_bufc    = 1.0/(1.0+Buf_C*Kbuf_C/(Y(3)+Kbuf_C)^2.0);
Ca_SR_bufSR = 1.0/(1.0+Buf_SR*Kbuf_SR/(Y(2)+Kbuf_SR)^2.0);

%% The Contractile Element (Calibrated Tran et al. 2017: https://doi.org/10.1113/JP274680)
% Model parameters

% Metabolite reference concentrations
Pi_ref = 2; %mM             Def: 2
MgATP_ref = 5; %mM          Def: 5
MgADP_ref = 0.036; %mM      Def: 0.036

% MgADP dissociation constant
kdADP = 0.004; % uM

% Set metabolite concentrations and pH level
MgATP = 5; %mM              Def: 5
MgADP = 36e-3; %mM          Def: 36e-3
Pii = 2; %mM                Def: 2
pH = 7.15; %                Def: 7.15                 

% Competitive H binding parameters
kdHCa = 1e-5;
m = 1;
pH_ref = 7.15;              %Def: 7.15
Href = 10^(-pH_ref)*1e3;

% Ischemic scenario
    if ischemiaFlag == 1  && myFlags(7)==1 % SEV1 (5 min) % Hypoxia on CE OFF:0 / ON:1
        MgATP = 4.2;
        MgADP = 0.33;
        Pii = 15.5;
        pH = 6.5;
    elseif ischemiaFlag == 2 && myFlags(7)==1   % SEV2 (10 min) % Hypoxia on CE OFF:0 / ON:1
        MgATP = 3.1;
        MgADP = 0.25;
        Pii = 22;
        pH = 6.3;
    end

% Gas constant
R = 8.314;
TmpC = 37;

% Sarcomere geometry (um)
SL_max = 2.4;
SL_min = 1.4;
L_thick = 1.65;
L_hbare = 0.1;
L_thin = 1.2;

% Temperature dependence
%TmpC = 25;
Q_kon = 1.5;
Q_koff = 1.3;
Q_kn_p = 1.6;
Q_kp_n = 1.6;
Q_fapp = 6.25;
Q_gapp = 2.5;
Q_hf = 6.25;
Q_hb = 6.25;
Q_gxb = 6.25;

% Species parameter
Xsp = 0.2; 

% XB density
rho = 0.25;                 % rho = 0.25;

% No of ATP consumed per XB 
phi = 1;

% Ca binding to troponin and thin filament
% fmyo = [kon_iso, koffL_iso, koffH_iso, perm50_iso, fapp_iso, gxb_iso, hb_iso, PExp_titin_iso, PCon_titin_iso];
kon = ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(1))*1.25*50e3 * GAC(23);   % (mM-1 s-1) 
koffL = ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(2)) * 200; % (s-1) Def: 250e-3
koffH = ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(3)) * 25;  % (s-1)
perm50 = ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(4)) * 0.6 * GAC(24); 
n_perm = 0.77*15;
kn_p = 550; % (s-1) % 
kp_n = 50; % (s-1)  % 
koffmod = 0.5;        % mod to change species

% Thin filament regulation and XB cycling
fapp = 1.0*500 * ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(5)) * GAC(25); % (s-1)
gapp = 1.0*70;  % (s-1)
gslmod = 1.0*6;
hf = 1*2000;  % (s-1)
hfmdc = 5;
hb = ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(7)) * 1.0*400;   % (s-1)
gxb = ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(6)) * 1.0*70;   % (s-1)
sigma_p = 8;
sigma_n = 1;   

% Mean strain of strongly-bond states
x0 = 0.007; % (um)
Psi = 2;

% Normalised active and passive force
SL_rest = 1.9;  % (um)
PCon_titin = 1.0*0.002; %(normalized force) Def: 0.002
PExp_titin = 1.0*10; %(um-1)
SL_collagen = 2.25; %(uM)
PCon_collagen = 1*0.02; %(normalized force)
PExp_collagen  = 1*70; %(um-1)

% Calculation of complete muscle response
mass = 4e-1*5e-5; % for Rat (normalised force s^2 um-1)
viscosity = 0.003;  % (normalised force s um-1)
    
    % Assign the state variables
    Nxb = Y(24);  % Non-permissive fraction in overlap region
    XBpreR = Y(25);   % Fraction of XBs in pre-rotated state (Stiffness generation)
    XBpostR = Y(26);  % Fraction of XBS in post-rotated state (force generation)
    x_XBpreR = Y(27);    % Strain of XBs in the pre-rotated state
    x_XBpostR = Y(28);  % Strain of XBs in the post-rotated state
    TropCaL = Y(29); % Ca bound to the low affinity troponin site
    TropCaH = Y(30); % Ca bound to the high affinity troponin site
    IntegF = Y(31); % Integral of force
    SL = Y(32); % Sarcomere length
 
    Afterload = 0; % Afterload for work-loop contraction. Def: 0
    loop = 0; % Boolean to indicate if it is a work-loop contraction
    preload_SL = 1.9; % Initial SL for work-loop contractions
    T_loop(1) = 2000; % Time for start of relaxation phasea (end of isotonic)
    T_loop(2) = 2000; % Tme for end of relaxation phase 
    passive = 1;      % Boolean for presence of passive force
    kxb = 13.1047;    % Def: 13.1047;

    % Prevent over extension
    if SL > preload_SL
        dSL = 0;
    end

  % *****************************************************    
    % Ca binding to troponin
    
    % Call function to get Ca
    Cai = Y(3); %mM
    
    % convert pH to mM
    H = 10^(-pH)*1e3; % H concentration in mM
    
    % Ca binding to troponin C
    konT = kon*1*Q_kon^((TmpC-37)/10);
    konT_app = konT*(kdHCa^m + Href^m)/(kdHCa^m + H^m);
    
    koffLT = koffL*koffmod*1*Q_koff^((TmpC-37)/10); 
    koffHT = koffH*koffmod*1*Q_koff^((TmpC-37)/10);
    
    d_TropCaL = konT_app*Cai*(1-TropCaL) - koffLT*TropCaL;
    d_TropCaH = konT_app*Cai*(1-TropCaH) - koffHT*TropCaH;
    
    
    % *****************************************************
    % Thin filament activation rates
    
    % Sarcomere geometry
    sovr_ze = min(L_thick/2, SL/2);
    sovr_cle = max(SL/2 - (SL-L_thin),L_hbare/2);
    L_sovr = sovr_ze - sovr_cle; % Length of single overlap region
    
    % Overlap fraction for thick filament
    sov_thick = L_sovr*2/(L_thick - L_hbare);
    % Overlap fraction for thin filament
    sov_thin = L_sovr/L_thin;
    
    TropReg = ((1-sov_thin)*TropCaL + sov_thin*TropCaH) * ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(12));
    permtot = (1/(1+(perm50/TropReg)^n_perm))^0.5;
    permtot_p_n = min(100,1/permtot);
    
    % Rate constants governing the transition btw Permissive and Non
    kn_pT = kn_p*permtot*Q_kn_p^((TmpC-37)/10);
    kp_nT = kp_n*permtot_p_n*Q_kp_n^((TmpC-37)/10);

    
    % *****************************************************
    % Cross-bridge cycling rates
    
    % Pxb to XBpreR
    fappT = fapp*Xsp*Q_fapp^((TmpC-37)/10);

    % Pi-dependent transition XBpreR to Pxb
    gappslmd = 1 + (1-sov_thick)*gslmod;
    gappT = gapp*gappslmd*Xsp*Q_gapp^((TmpC-37)/10);
    gappT_true = gappT/Pi_ref;  % True first order rate constant
    
    % XBpreR to XBpostR
    hfmd = exp(-sign(x_XBpreR)*hfmdc*(x_XBpreR/x0)^2);
    hfT = hf*hfmd*Xsp*Q_hf^((TmpC-37)/10);
    
    % H-dependent transition XBpostR to XBpreR
    hbT = hb*Xsp*Q_hb^((TmpC-37)/10);
    hbT_true = hbT/Href;
    hbT_app = ((kdADP+MgADP_ref)/MgADP_ref)*(MgADP/(kdADP+MgADP))*hbT_true;  

    % MgATP-dependent transition from XBpostR to Pxb
    if (x_XBpostR < x0)
        gxbmd = exp(sigma_p*((x0-x_XBpostR)/x0)^2);
    else
        gxbmd = exp(sigma_n*((x0-x_XBpostR)/x0)^2);
    end
    
    gxbT = gxb*max(gxbmd,1)*Xsp*Q_gxb^((TmpC-37)/10);
    gxbT_true = gxbT/MgATP_ref;
    gxbT_app = ((kdADP+MgADP_ref)/(kdADP+MgADP))*gxbT_true;

    % Pxb to XBpostR - Introduced for thermodynamic efficiency
    G0 = -29600 + log(10)*R*(TmpC+273)*(-log10(1e-7));
    K_MgATP = exp(-G0./(R.*(TmpC+273)))*1e6;  % 1e6 to convert from M2 to mM2
    fxbT = (kdADP*fappT*hfT*gxbT_true)/(gappT_true*hbT_true*K_MgATP);
    fxb = (kdADP*fapp*hf*(gxb/MgATP_ref))/((gapp/Pi_ref)*(hb/Href)*K_MgATP); % Used for calculating max
    
    ap1 = fappT;
    ap2 = hfT * ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(13));
    ap3 = MgATP*gxbT_app * ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(14));
    am1 = Pii*gappT_true;
    am2 = H*hbT_app * ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(15));
    am3 = fxbT;
    F1 = 1 * ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(10));
    F2 = 1 * ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(11));
    
    % *****************************************************
    %  Cross_bridge transitions
    
    % Update all RUs   
    Pxb = 1 - Nxb - XBpreR - XBpostR;
    d_Nxb = -kn_pT*Nxb + kp_nT*Pxb;
    
    d_XBpreR = ap1*F1*Pxb - am1*F2*XBpreR - ap2*XBpreR + am2*XBpostR;
    d_XBpostR = ap2*XBpreR + am3*Pxb - am2*XBpostR - ap3*XBpostR;
    
    ATPase = rho*phi*sov_thick*(ap3*XBpostR - am3*Pxb);
    
    % *****************************************************
    % Mean strain of strongly-bound states
           
    % Steady State occupancies of the bound states - Duty fractions
    Sum_duty = am3*am2 + ap3*ap1 + am2*ap1 + ap1*ap2 + am3*am1 + am3*ap2...
        + ap2*ap3 + am3*am1 + ap3*am1;

    dutyPreR = (am3*am2 + ap3*ap1 + am2*ap1)/Sum_duty;
    dutyPostR = (ap1*ap2 + am3*am1 + am3*ap2)/Sum_duty;
    
    % No shortening during isotonic period
    % This is to hold the SL at max contraction to get a loop
    if (time>T_loop(1)) & (time<T_loop(2)) & loop
        dSL = 0;
    end

    % Compute the active force
    F_active = sov_thick*(x_XBpreR*XBpreR + x_XBpostR*XBpostR);

    % Maximal state occupancies under optimal conditions
    max_XBpreR = (fapp*(hb+gxb)+fxb*fapp)/(gxb*hf + fapp*hf + gapp*hb + ...
        gapp*gxb + fapp*hb + fapp*gxb + fxb*(hb+gapp+hb));

    max_XBpostR = (fapp*hf+fxb*(gapp+hb))/(gxb*hf + fapp*hf + gapp*hb + ...
        gapp*gxb + fapp*hb + fapp*gxb + fxb*(hb+gapp+hb));

    % Factor to normalise force
    FnormD = x0*max_XBpostR;
    
    % Normalised Active force
    Fnorm = F_active/FnormD;
    
    % Normalised Passive force
    F_passive = passiveForces(SL,time,tDrugApplication,fmyo);
    F_after = Afterload;
    
    % No shortening during isotonic period
    % This is to hold the SL at max contraction to get a loop
    if (time>T_loop(1)) & (time<T_loop(2)) & loop
        dSL = 0;
    end
    
    % Difference in force
    d_Force =  F_after -  Fnorm - F_passive;
    
    % Total normalised force
    F_total = Fnorm + abs(F_passive);
    
    % Used in Fevents script 
    F_net = F_total - F_after;
    
    % SL dynamics
    
    dSL = contrflag*((IntegF + viscosity*(SL_rest - SL))/mass);
    %F_passive = passiveForces(SL,passive);
    %dIntegF = (-F_passive-Fnorm+Afterload);
    % Rate of change of the average distortions
    d_x_XBpreR = dSL/2 + (Psi/dutyPreR)*(-ap1*x_XBpreR) + am2*(x_XBpostR-x0-x_XBpreR);
    d_x_XBpostR = dSL/2 + (Psi/dutyPostR)*(ap2*(x_XBpreR + x0 - x_XBpostR));
    
    % Ca buffering by low-affinity troponin C (LTRPNCa)
    Trop_conc = 70e-3;       % (mM) troponin concentration
    
    FrSBXB    = (XBpostR+XBpreR)/(max_XBpostR + max_XBpreR);         %Done
    dFrSBXB   = (d_XBpostR +d_x_XBpreR)/(max_XBpostR + max_XBpreR);  %Done

    dsovr_ze  = -dSL/2*heav(L_thick-SL);
    dsovr_cle = -dSL/2*heav((2*L_thin-SL)-L_hbare);
    dlen_sovr = dsovr_ze-dsovr_cle;
    dSOVFThin = dlen_sovr/L_thin;

    dTropTot= Trop_conc*(-dSOVFThin*TropCaL+(1-sov_thin)*d_TropCaL + ...        
    dSOVFThin*(FrSBXB*TropCaL+(1-FrSBXB)*TropCaL) + ...
    sov_thin*(dFrSBXB*TropCaL+FrSBXB*d_TropCaH-dFrSBXB*TropCaL+...
    (1-FrSBXB)*d_TropCaL));
     
    % Assembling the derivative vector
    dY(24:32,1) = [d_Nxb, d_XBpreR, d_XBpostR, d_x_XBpreR, d_x_XBpostR, d_TropCaL, d_TropCaH, d_Force, dSL];
    
    JltrpnLH = dTropTot;
    JCaBMyo = JltrpnLH;
    frc = kxb*Fnorm;           % 2024 Simulations (only reporting active)
    frc = kxb*F_total;         % Reporting active + passive
    Lsarc = SL;
    cvelo   = dSL;

    if contrflag == 2 % No contraction
            JltrpnLH = 0;
            JCaBMyo = 0;
            frc = 0;
            Lsarc = SL_rest;
            cvelo   = 0;
            ATPase = 0;
    end


%% Oxygen dynamics

OBL = 0.1772;            % Source/bath O2 concentration in normal condition (physoxia) 0.133 mM (https://doi.org/10.1074/jbc.M116.751826) | OBL = 0.17956 for 17% or 0.1625 for 15%
alp = 1.1626e-01;       % conversion factor Cm/(F*Vc*1e-18): A/F to mM/s Paci et al., 2018 and Wei et al., 2014
bta = (1-z)*OBL;        % Physoxia = 0, SEV1 = 0.5855 representing 20% INaK & IpCa inhibition, SEV2 = 0.7174 representing 31% INaK & IpCa inhibition. bta = 1-O2.

D_m2s   = 2.18e-05 * 1e-4;       % 2.18e-05 cm^2/s -> 2.18e-09 m^2/s
H_eff_m = 6.22e-3;               % the measured whole medium height (6.22 mm)
%H_eff_m = 1.0e-4;               % effective control volume height (100 um)
%eps     = D_m2s / (H_eff_m^2);  % (s^-1) diffusion rate based on Fick's law ≈ 5.63e-05 s^-1
eps = 2.56e-04;                  


ATP_basal = 0.01;          % Non-ATP-linked mitochondrial respiration (proton leak/maintenance). ATP_basal = 0.01;
J_SERCA = i_up;

O2s = OBL-bta*(0.5*(1+tanh((0.06)*(time-300))))+bta*(0.5*(1+tanh((0.06)*(time-900)))); 

% Stoichiometric O2-per-ATP factor
phiATP = 0.20;   % O2 (mM) per ATP (mM)

% --- Convert processes to ATP rates (mM/s ATP) ---
% alp converts A/F -> mM/s of monovalent charge. For NaK, 1 ATP per net +1 charge.
ATP_NaK   = alp * abs(i_NaK) * 1.0;      % 1 ATP per pump cycle ≈ 1 net charge

J_PMCA_Ca = (alp/2) * abs(i_PCa);        % i_PCa: Ca2+ current (A/F) -> Ca flux: use z=2mM/s Ca
ATP_PMCA  = J_PMCA_Ca * 1.0;             % 1 ATP per Ca2+ extruded
ATP_SERCA = 0.5 * J_SERCA;               % 1 ATP per 2 Ca2+ pumped (J_SERCA in mM/s)
ATP_XB    = ATPase;                      % the myosin ATPase in mM/s ATP

% Optional basal/maintenance (ATP not tied to specific fluxes)
ATP_leak  = ATP_basal;                   % e.g., small non-ATP-linked mitochondrial respiration (proton leak/maintenance).

% 1) Raw OCR terms (positive mM/s O2)
OCR_XB_raw    = phiATP .* ATP_XB;        
OCR_SERCA_raw = phiATP .* ATP_SERCA;
OCR_NaK_raw   = phiATP .* ATP_NaK;
OCR_PMCA_raw  = phiATP .* ATP_PMCA;

% Leak: since ATP_basal is NON-ATP-linked proton leak, we treat it as O2 directly
OCR_leak_raw  = ATP_basal;               % mM/s O2  

OCR_nonXB_raw = OCR_SERCA_raw + OCR_NaK_raw + OCR_PMCA_raw + OCR_leak_raw;

alpha = 0.008;                         

% Apply alpha to NON-XB pathways (preserves their relative proportions)
OCR_XB    = OCR_XB_raw;
OCR_SERCA = alpha .* OCR_SERCA_raw;
OCR_NaK   = alpha .* OCR_NaK_raw;
OCR_PMCA  = alpha .* OCR_PMCA_raw;
OCR_leak  = alpha .* OCR_leak_raw;


OCR_XB    = gma .* OCR_XB;
OCR_SERCA = gma .* OCR_SERCA;
OCR_NaK   = gma .* OCR_NaK;
OCR_PMCA  = gma .* OCR_PMCA;
OCR_leak  = gma .* OCR_leak;

OCR_nonXB = OCR_NaK + OCR_PMCA + OCR_SERCA + OCR_leak;
% Finish state equation
OCR_terms = OCR_XB + OCR_SERCA + OCR_NaK + OCR_PMCA + OCR_leak;   % +ve consumption
OCR       = -OCR_terms;


% ---------- mapping per-cell OCR to pericellular medium ----------
Awell = 3.217e-05;                                % well area (m^2)
delta_m   = 1.0e-4;                               % m, represents a thin pericellular film (e.g., 50–200 µm)
Vc_L      = 8800 * 1e-15;                         % cell volume [L] (8800 µm^3)
V_peri_L  = 200e-6;                               % entire well volume [L] (200 µL from Emilia's table)
%V_peri_L  = Awell * delta_m * 1000;              % pericellular control volume m^3 -> L
N         = 5.0e4;                                % cells/well (use Emilias's experiment's 50,000)
scale_peri = (N * Vc_L) / V_peri_L;               % dimensionless ~0.0022.


dY(33,1)  = OCR * scale_peri + eps*(O2s - Y(33));

O2O = max(Y(33), 0);              
do2dt = dY(33,1);
PRF =  eps*(O2s-Y(33));

% Calculating pmol/min readouts to enable validation with in vitro data
r_mMs_cell = -OCR;                         % mM/s per cell (positive = consumption)

% Cell volume used in the model:
Vc_um3    = 8800.0;                    % um^3
V_cell_mL = Vc_um3 * 1e-12;            % 1 um^3 = 1e-12 mL  -> 8.8e-9 mL
K_cell    = V_cell_mL * 60e6;          % = 0.528 pmol/min per (mM/s)


% Instantaneous amount rates
pmol_min_cell = r_mMs_cell .* K_cell;        % pmol/min per cell
pmol_min_well = pmol_min_cell .* N;          % pmol/min per well


%% Ionic concentrations
%Nai
dY(18, 1)   = -Cm*(i_Na+i_NaL+i_b_Na+3.0*i_NaK+3.0*i_NaCa+i_fNa)/(F*Vc*1.0e-18);
%Cai
%dY(3, 1)    = Cai_bufc*(i_leak-i_up+i_rel-(i_CaL+i_b_Ca+i_PCa-2.0*i_NaCa)*Cm/(2.0*Vc*F*1.0e-18)-JCaBMyo);
dY(3,1)      = Cai_bufc * (i_leak - i_up + i_rel - ((i_CaL + i_b_Ca + i_PCa - 2.0*i_NaCa) * Cm / (2.0 * Vc * F * 1.0e-18)) - JCaBMyo);
%caSR
dY(2, 1)    = Ca_SR_bufSR*Vc/V_SR*(i_up-(i_rel+i_leak));

%% Stimulation
i_stim_Amplitude 		= 5.5e-10;   % 7.5e-10 or 5.5e-10  % ampere (in stim_mode)
i_stim_End 				= 20000.0;   % second (in stim_mode)    
i_stim_PulseDuration	= 0.005;     % second (in stim_mode)
i_stim_Start 			= 0.0;       % second (in stim_mode)
if classicVSoptic == 0
    i_stim_frequency 	= 60.0;   % 
elseif classicVSoptic == 1
    i_stim_frequency 	= 30.0;   % 
    %i_stim_frequency 	= 60.0;   % 
end
stim_flag 				= stimFlag;   % dimensionless (in stim_mode)
%i_stim_Period 			= 60.0/i_stim_frequency;
i_stim_Period = 1/freq;

if stimFlag>1
error('hiMCES: stimFlag must be only 0 (spontaneous) or 1 (paced)');
end

if ((time >= i_stim_Start) && (time <= i_stim_End) && (time-i_stim_Start-floor((time-i_stim_Start)/i_stim_Period)*i_stim_Period <= i_stim_PulseDuration))
    i_stim = stim_flag*i_stim_Amplitude/Cm;
else
    i_stim = 0.0;
end

%% Membrane potential
I_bias = 0;

dY(1, 1) = -(i_K1+i_to+IKatp+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_NaL+i_NaCa+i_PCa+i_f+i_b_Na+i_b_Ca-i_stim+I_bias);

%% Output variables
IK1     = i_K1;
Ito     = i_to;
IKr     = i_Kr;
IKs     = i_Ks;
ICaL    = i_CaL;
INaK    = i_NaK;
INa     = i_Na;
INaCa   = i_NaCa;
IpCa    = i_PCa;
If      = i_f;
IbNa    = i_b_Na;
IbCa    = i_b_Ca;
Irel    = i_rel;
Iup     = i_up;
Ileak   = i_leak;
Istim   = i_stim;
INaL    = i_NaL;
JCaB    = JCaBMyo;
psph    = phosph;
vcy     = v_cycle;

dati = [INa, If, ICaL, Ito, IKs, IKr, IK1, INaCa, INaK, IpCa, IbNa, IbCa, Irel, Iup, Ileak, Istim, ...
    E_K, E_Na, INaL, F_passive, frc, ATPase, Lsarc, cvelo, JCaB, psph, IKatp, ronak, O2s, O2O, vcy, do2dt, ...
    finhib, Ko, ro1, ro2, fkatp, kon, OCR, OCR_XB, OCR_nonXB, PRF, pmol_min_cell, ...
    pmol_min_well, OCR_NaK, OCR_leak, OCR_SERCA, OCR_PMCA];
%===============================================================================
% End
%===============================================================================
