function [ n ] = getRefractiveIndex( wavelength, varargin )
%getRefractiveIndex(L,S,T) Returns the refractive index of sea water given
%the wavelength (in nanometers), salinity (in parts per thousand), and
%temperature (in degrees Celsius).
%
%getRefractiveIndex(L,S) Returns the refractive index of water assuming a
%temperature of 19 degrees Celsius.
%
%getRefractiveIndex(L) Returns the refractive index assuming a salinity of
%35 parts per thousand and a temperature of 19 degrees Celsius.
%
% Robert Foster 2014-10-08

%Default values if no inputs are given
defaultSalinity     = 35;   %Parts per thousand
defaultTemperature  = 19;   %Degrees Celsius

%==========================================================================
%%%%%%%%%%%%%%%%%%%%%% C H E C K  I N P U T S %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================

%Check to make sure the wavelength input is a vector
if(~isvector(wavelength))
    error('Wavelength input must be a vector! (Size must be [1 n] or [n 1]!)\n');
end

%If only wavelength is given, assume default salinity and temperature
if (nargin == 1)
    salinity = defaultSalinity.*ones(size(wavelength));
    temperature_C = defaultTemperature.*ones(size(wavelength));
    
    %If wavelength and salinity are given, assume default temperature
elseif (nargin == 2)
    temperature_C = defaultTemperature.*ones(size(wavelength));
    
    %If passed a scalar value, make it a vector the same size as wavelength
    salinity = varargin{1}.*ones(size(wavelength));
    
    %When all arguments are given
elseif(nargin >= 3)
    %If passed a scalar value, make it a vector the same size as wavelength
    salinity = varargin{1}.*ones(size(wavelength));
    temperature_C = varargin{2}.*ones(size(wavelength));
end
%==========================================================================
%%%%%%%%%%% C A L C U L A T E  R E F R A C T I V E  I N D E X %%%%%%%%%%%%%
%==========================================================================

%Equation from Quan & Fry, 1995 "Empirical equation for the index 
%                               of refraction of seawater", Applied Optics

%Constants
n0 = 1.31405;
n1 = 1.779e-4;
n2 = -1.05e-6;
n3 = 1.6e-8;
n4 = -2.02e-6;
n5 = 15.868;
n6 = 0.01155;
n7 = -0.00423;
n8 = -4382;
n9 = 1.1455e6;

%Calculate refractive index 
n = n0 + ...
    salinity.*(n1 + n2.*temperature_C + n3.*(temperature_C.^2) ) + ...
    n4.*(temperature_C.^2) + ...
    (n5 + n6.*salinity + n7.*temperature_C)./wavelength + ...
    n8./(wavelength.^2) + ...
    n9./(wavelength.^3);

end

