function [DepolarizationFactor] = getAtmosphericDepolarizationFactor(wav)
%Depolarization factor of the atmosphere is wavelength dependent
%Table of values from Bucholtz (1995)
%Wavelength ranges from 200 nm to 1000 nm

%Input wavelength can be either in microns or nanometers

if (wav <= 1000 && wav >= 200)
    wav = wav/1000;
end

if (wav < 0.2 || wav > 1)
    error('Wavelength must be between 200 and 1000 nanometers (0.2 and 1 microns)')
end

adf_wav = [0.2,0.205,0.21,0.215,0.22,0.225,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,...
            0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,...
                0.8,0.85,0.9,0.95,1];
adf = [4.545,4.384,4.221,4.113,4.004,3.895,3.785,3.675,3.565,3.455,3.4,3.289,3.233,3.178,...
        3.178,3.122,3.066,3.066,3.01,3.01,3.01,2.955,2.955,2.955,2.899,2.842,2.842,2.786,...
         2.786,2.786,2.786,2.73,2.73,2.73,2.73,2.73];

adf = adf./100;

DepolarizationFactor = interp1(adf_wav, adf, wav, 'linear');
end