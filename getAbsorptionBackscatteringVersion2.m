function [a_phy,a_CDOM,a_NAP,b_phy,bb_phy,b_NAP,bb_NAP,a_water,bb_water,b_water] = getAbsorptionBackscatteringVersion2(CHL,wav)
%% Absorption and Backscattering for Case 1 Waters

% Based off Gilersons method in 
% Bio-optical Modeling of Sun-Induced Chlorophyll-a Fluorescence, Ch 7, Gilerson-Huot 

% As well as parts from Light and Water, Ch 3, Mobley

% Backscatter ratio of pure water taken from www.oceanopticsbook.info

% Uses Spline interpolation

%% Reading Files
wav_pico = xlsread('PicoandMicroAbsorption.xlsx','C1:C201');
a_pico_file = xlsread('PicoandMicroAbsorption.xlsx','B1:B201');
a_micro_file = xlsread('PicoandMicroAbsorption.xlsx','A1:A201');

wav_water = xlsread('PureWaterIOPs.xls','A2:A162');
a_water_file = xlsread('PureWaterIOPs.xls','B2:B162');
bb_water_file = xlsread('PureWaterIOPs.xls','C2:C162');

%% Interpolating Data
a_pico = interp1(wav_pico,a_pico_file,wav,'spline','extrap');
a_micro = interp1(wav_pico,a_micro_file,wav,'spline','extrap');

a_water = interp1(wav_water,a_water_file,wav,'spline','extrap');
bb_water = interp1(wav_water,bb_water_file,wav,'spline','extrap');

%% Water Scattering

B = 0.5;        %Backscattering Ratio - bb/b

b_water = bb_water./B;

%% Phytoplankton Absorption
a_star_phy_lam = 0.042*CHL^(-0.2);
a_pico_lam = interp1(wav_pico,a_pico_file,443);
a_micro_lam = interp1(wav_pico,a_micro_file,443);

Sf = (a_star_phy_lam - a_micro_lam)/(a_pico_lam - a_micro_lam);
a_star_phy = Sf.*a_pico + (1-Sf).*a_micro;
a_phy = a_star_phy.*CHL;

a_phy_443 = interp1(wav,a_phy,443);

%% CDOM Absorption

x = -0.014;     %Range from -0.014 to -0.019 (-0.014 typical)
a_CDOM_443 = 0.7*a_phy_443;

a_CDOM = a_CDOM_443.*exp(x.*(wav-443));

%% NAP Absorption

y = -0.011;     %Range from -0.006 to -0.014 (-0.011 typical)
a_NAP_443 = 0.56*a_phy_443;

a_NAP = a_NAP_443.*exp(y.*(wav-443));

NAP = a_NAP_443/0.035;

%% Phytoplankton Scattering and Backscattering

c_phy_550 = 0.3*CHL^0.57;
c_phy = c_phy_550.*(550./wav).^0.8;

b_phy = c_phy - a_phy;
bb_phy =  0.006.*b_phy;

%% NAP Scattering and Backscattering

b_NAP_550 = 0.75*NAP;
b_NAP = b_NAP_550.*(550./wav).^0.8;
bb_NAP = 0.0183.*b_NAP;

end