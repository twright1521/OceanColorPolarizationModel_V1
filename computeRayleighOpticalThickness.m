function [ RayleighOpticalThickness, ...
           RayleighVolScatCoeff,     ...
           RayleighCrossSection ]    ...
               = computeRayleighOpticalThickness( wavelength_nm, varargin )
% computeRayleighOpticalThickness 
%  Calculates the Rayleigh optical thickness, Rayleigh Volume Scattering
%  Coefficient, and the Rayleigh Scattering Cross Section.  Based on the
%  analytic fitted functions in:
%       'Rayleigh-scattering calculations for the terrestrial atmosphere'
%       Anthony Bucholtz, Applied Optics, 1995
%
% Required Inputs:
%    (1) Wavelength       - [nanometers].  This function requires the
%                           wavelength of light (in vacuum).
%
% Optional Inputs:
%    (2) Atmosphere Model - 'Tropical', 'MidLatSummer', 'MidLatWinter', 
%                           'SubArcticSummer', 'SubArcticWinter',
%                           '1962Standard'.  Assumes '1962Standard' if not
%                           specified.
%    (3) Pressure         - [mbars].  Pressure where you want the
%                           calculation to take place.  Assumes 1013.25 if
%                           not specified.
%    (4) Temperature      - [Kelvin]. Temperature where you want the
%                           calculation to take place.  Assumes 288.15 K if
%                           not specified.
%
%   The function can accept an array of wavelengths, but only one
%   atmosphere, pressure, and temperature input which is applied to all
%   wavelengths.  This may change in the future if I get ambitious.
%
%   Robert Foster
%   rfoster01@citymail.cuny.edu
%   2015-03-17
%


%Reference Pressure and Temperature
P_ref = 1013.25; %mbars
T_ref = 288.15;  %Kelvin

%Initialize variables
[A, B, C, D] =  deal(nan(size(wavelength_nm)));

%% 
%==========================================================================
%%%%%%%%  C O M P U T E  R A Y L E I G H  C R O S S  S E C T I O N %%%%%%%%
%==========================================================================
%Parameters
wavelength_um   = wavelength_nm ./ 1000;
wav_lt500       = wavelength_um <= 0.5;
wav_gt500       = wavelength_um > 0.5;

%Check to see if any wavelengths are outside the valid range 
if(any(wavelength_um < 0.2))
   error('This function does not work for wavelengths less than 200nm!\n');
end

%Fitting Coefficients  (Table 3)
%For wavelengths < 0.5um
A(wav_lt500) = 3.01577e-28;
B(wav_lt500) = 3.55212;
C(wav_lt500) = 1.35579;
D(wav_lt500) = 0.11563;

%For wavelengths > 0.5um (Table 3)
A(wav_gt500) = 4.01061e-28;
B(wav_gt500) = 3.99668;
C(wav_gt500) = 1.10298e-3;
D(wav_gt500) = 2.71393e-2;

%Compute Rayleigh Cross sections
RayleighCrossSection = ...
                       A.*wavelength_um.^ ...
                      (-1.*( B ....
                           + C.*wavelength_um ...
                           + D./wavelength_um ...
                           ));

%%                                 
%==========================================================================
%%%%%%% R A Y L E I G H  V O L U M E  S C A T T E R I N G  C O E F F %%%%%%
%==========================================================================

%Evaluate Inputs
%No inputs other than wavelength, use standard atm, temp and pressure
if(nargin == 1)
    Atm     = '';
    P_meas  = P_ref;
    T_meas  = T_ref;
    warning('%s Using P=%6.1fmbar and T=%6.1fK.\n', ...
            'No Pressure or Temperature specified.', P_meas, T_meas);
        
%Only Atmosphere is given
elseif(nargin == 2)
    Atm     = varargin{1};
    P_meas  = P_ref;
    T_meas  = T_ref;
    warning('%s Using P=%6.1fmbar and T=%6.1fK.\n', ...
            'No Pressure or Temperature specified.', P_meas, T_meas);
        
%Atmosphere and Pressure are given        
elseif(nargin == 3)
    Atm     = varargin{1};
    P_meas  = varargin{2};
    T_meas  = T_ref;
    warning('No Temperature specified. Using T=%6.1fK.\n', T_meas);
    
%Atmosphere Pressure and temperature are given    
elseif(nargin == 4)
    Atm     = varargin{1};
    P_meas  = varargin{2};
    T_meas  = varargin{3};
end

%For wavelengths < 0.5um (Table 3)
A(wav_lt500) = 7.68246e-4;

%For wavelengths > 0.5um (Table 3)
A(wav_gt500) = 10.21675e-4;

%Coefficients B, C, and D do not change.
%Compute Rayleigh Volume Scattering Coefficients
RayleighVolScatCoeff = ...
                       A.*wavelength_um.^ ...
                       (-1.*( B ....
                            + C.*wavelength_um ...
                            + D./wavelength_um ...
                            ));
                       
%Adjust the Volume Scattering Coeff for pressure and temperature. (Eqn. 10)
RayleighVolScatCoeff =  RayleighVolScatCoeff .* ...
                       (P_meas./P_ref) .* ...
                       (T_ref./T_meas);
                   
%%                                 
%==========================================================================
%%%%%%%%%%%% R A Y L E I G H  O P T I C A L  T H I C K N E S S  %%%%%%%%%%%
%==========================================================================                   

%Get the atmospheric model surface pressure (mbar) temperature (K), and A
%fitting coefficient (Table 5)
switch Atm
    case 'Tropical'
        %Tsurf_model     = 300;
        Psurf_model     = 1013;
        A(wav_lt500)    = 6.52965e-3;
        A(wav_gt500)    = 8.68094e-3;       
    case 'MidLatSummer'
        %Tsurf_model     = 294;
        Psurf_model     = 1013;
        A(wav_lt500)    = 6.51949e-3;
        A(wav_gt500)    = 8.66735e-3;
    case 'MidLatWinter' 
        %Tsurf_model     = 272.2;
        Psurf_model     = 1018;
        A(wav_lt500)    = 6.53602e-3;
        A(wav_gt500)    = 8.68941e-3;
    case 'SubArcticSummer'
        %Tsurf_model     = 287;
        Psurf_model     = 1010;
        A(wav_lt500)    = 6.48153e-3;
        A(wav_gt500)    = 8.61695e-3;
    case 'SubArcticWinter'
        %Tsurf_model     = 257.1;
        Psurf_model     = 1013;
        A(wav_lt500)    = 6.49997e-3;
        A(wav_gt500)    = 8.64145e-3;
    case '1962Standard' 
        %Tsurf_model     = 288.1;
        Psurf_model     = 1013;
        A(wav_lt500)    = 6.50362e-3;
        A(wav_gt500)    = 8.64627e-3;
    otherwise
        warning(...
            '%s %s\nValid options are %s, %s, %s, %s, %s, and %s.\n', ...
            'No Atm. Model specified, or wrong spelling. ', ...
            'Using standard 1962 Model.',   ...
            'Tropical',                     ...
            'MidLatSummer',                 ...
            'MidLatWinter',                 ... 
            'SubArcticSummer',              ...
            'SubArcticWinter',              ...
            '1962Standard');
        %Tsurf_model     = 288.1;
        Psurf_model     = 1013;
        A(wav_lt500)    = 6.50362e-3;
        A(wav_gt500)    = 8.64627e-3;
end

%Coefficients B, C, and D do not change.
%Compute Rayleigh Optical Thickness
RayleighOpticalThickness = ...
                       A.*wavelength_um.^ ...
                       (-1.*( B ....
                            + C.*wavelength_um ...
                            + D./wavelength_um ...
                            ));
                       
%Correct the Rayleigh Optical Thickness for pressure  (Eqn. 17)
RayleighOpticalThickness = ...
                    RayleighOpticalThickness .* (P_meas./Psurf_model);

