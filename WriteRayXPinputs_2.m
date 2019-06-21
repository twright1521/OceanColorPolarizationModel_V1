%% |||||||||||||||DISCRIPTION |||||||||||||||||||||||||||||||||||||||||||||
%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% This program creates Ray XP input files for designated wavelengths 
% and designated chlorophyl concentration.

% This program is specifically for a single chl concentration and multiple
% wind speeds.

% It also plots the absorption, scattering, and backscattering coefficients
% at the specified wavelengths and chlorophyl concentration
% if so desired.

clearvars

%% |||||||||||||||INFO TO CHANGE AT THE START OF EVERY SIMULATION |||||||||
%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sim_num = 4;        % Denotes the simulation folder number
                    % Determines which Simulation folder
                    % the files will be saved in
                    % Change this number for a new group of simulations

Version_num = 1;    % Used to denote different versions
                    % for the same chlorophyl and wind speed inputs.
                    % Best used when changing other inputs
                    % such as Sun and Receiver angles or other parameters.
                 
CHL = 30;        % Chlorophyl concentration in milligrams/m^3
                
% Specify the wavelengths to run or a range of wavelengths
wav = linspace(400,900,101); %[200 to 900] nanometers

[~,m] = size(wav);

% Would you like to plot absorption and scattering?

plotTag = false;        % Set to true if plotting is desired
                        % Set to false if undesired
%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% Wind speed for waves
WS = [0,3,5,10];       %Windspeed in m/s = 0 - no Windy waves [2.5 .. 15] Windy waves
                       
% Parameters to be varied for more dynamic results
SSA = 0.905;     %Single Scatter Albedo for near ocean atmosphere, [0.85 - 0.95] - Mean is 0.905 
InclDirectSolarBeam = 0; %Use Values greater than 0 to include direct light
Atm_Temp = 277;      %Atmospheric Temperature - Kelvin

% Astronomical Angles - Integers only
Sun_Pol = 30;       % Sun polar angle, degrees                   
Rec_Pol = 60;       % Receiver polar angle, degrees (integer)
Rec_Az = 90;        % Receiver azimuth relative to sun, degrees (integer) 

% Angle Files 
% Angles specified above must be in these files, respectively
sun_pol_file = '_Ang_30_45_60.amu';
rec_pol_file = '_PolarScale.amu';
rec_azm_file = '_AzimuthScale_0_90_180_270.azm';

% Constants 
depth_assump = 1000;    %Assumed depth
Delta_hsol = 0.039;     %Constant

%% ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  P R O G R A M  S T A R T %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fignum = 1;


for zeta = 1:length(WS)
    % Location to save input files
    sim_fol_name = sprintf('Simulation_%u',Sim_num);

    if (exist(sim_fol_name, 'dir') == 0)
        mkdir(sim_fol_name)
    end

    inp_fol_name = strcat('RayXP_CHL_',num2str(CHL),'_WS_',num2str(WS(zeta)),'_V_',num2str(Version_num));

    mkdir(strcat(sim_fol_name,filesep,inp_fol_name))
    mkdir(strcat(sim_fol_name,filesep,inp_fol_name,filesep,'Input Files'))
    mkdir(strcat(sim_fol_name,filesep,inp_fol_name,filesep,'Output Files'))

    %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    % Absorption and Scattering

    [a_phy,a_CDOM,a_NAP,b_phy,bb_phy,b_NAP,bb_NAP,a_water,bb_water,b_water] = getAbsorptionBackscatteringVersion2(CHL,wav);

    % For plotting purposes only
    if (plotTag == true)
        A_total = a_phy + a_CDOM + a_NAP + a_water;
        B_total = b_phy + b_NAP + b_water;
        Bb_total = bb_phy + bb_NAP + bb_water;
    end

    % Optical Properties of Hydrosol

    C_sol = a_phy + a_NAP + b_phy + b_NAP;
    SSA_hsol = (b_phy + b_NAP)./C_sol;

    TauMol_hsol = b_water * depth_assump;
    TauAbs_hsol = (a_water + a_CDOM) * depth_assump;
    TauSol_hsol = C_sol * depth_assump;

    % Wind Check

    if (WS(zeta) == 0)
        waveName = 'NULL';
        windTag = false;
    else 
        waveName = 'ISOTROPIC';
        windTag = true;
    end

    %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    % Writing the RayXP File

    for x = 1:m                                

        %HSL file
        if (wav(x) <= 477)
            HSL_str = 'PTH_mmVF.440';
        else
            if (wav(x) > 477 && wav(x) <= 532)
                HSL_str = 'PTH_mmVF.514';
            else 
                if (wav(x) > 532 && wav(x) <= 613)
                    HSL_str = 'PTH_mmVF.550';
                else 
                    if (wav(x) > 613)
                        HSL_str = 'PTH_mmVF.675';
                    end
                end
            end
        end
        RefractiveIndex = getRefractiveIndex(wav(x), 35, 19);      %RefractiveIndex
        %Create unique input filenames based on varying parameters
        filename = sprintf('Test_Script_%3u', wav(x));

        %Open file and write the RayXP header
        filepath = strcat(sim_fol_name,filesep,inp_fol_name,filesep,'Input Files',filesep,filename, '.in3');
        IFILE = fopen(filepath, 'w');

        %Write header information    
        fprintf(IFILE, 'IN3File\n');
        fprintf(IFILE, 'Wavelength: %f\n', wav(x)/1000); %Wavelength in micrometers
        fprintf(IFILE, 'RefractiveIndex: %f\n', RefractiveIndex);
        fprintf(IFILE, 'START_LAYERS:\n');

        %Values for AOT for TAUSOL
        %These values create a specific atmosphere taken from literature 
        %(exactly where? I dunno)
        aot = [0.28725,0.2305,0.1415,0.09625,0.055]; %aerosol optical thickness
        aotwave = [380,500,675,870,1020];

        %Get the Rayleigh Optical Thickness
        [TauRay, ~, ~] = computeRayleighOpticalThickness(wav(x), 'MidLatSummer',1013.25,Atm_Temp);

        %%%%%%%%%%%%%%%%%%%%%%%%%% SOURCE LAYER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SRCtype = 'LOWTRAN7'; % SRCtype:[NULL,LOWTRAN7]
                                    % NULL - ISun - ETSI, W/(m*m*micrometer)
                                    % LOWTRAN7 - ISun ignored, ETSI evaluates to LOWTRAN7 model data for assigned "Wavelength" 
        ISun = 0;                   % ISun >= 1.E-3, W/(m*m*micrometer)
        MoonRefl = 0;         % MoonRefl = [1.E-8..1] - Moon reflection.
                                    %If 1.E-8<=MoonRefl<=1, then EXTRA-TERRESTRIAL IRRADIANCE = MoonRefl*ETSI; 
        Par2 = 0;             % Par[2] - reserved
        Par3 = 0;             % Par[3] - reserved
        Par4 = 0;             % Par[4] - reserved
        Par5 = 0;             % Par[5] - reserved
        fprintf(IFILE, '%-15s%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n', '//SOURCE:', 'SRCtype', 'ISun', 'MoonRefl', 'NotUsed', 'NotUsed', 'NotUsed', 'NotUsed');
        fprintf(IFILE, '%-15s%-10s%-9.7g %-9.7g %-9.7g %-9.7g %-9.7g %-9.7g \n\n', '  SOURCE:', SRCtype, ISun, MoonRefl, Par4, Par3, Par4, Par5);    

        %%%%%%%%%%%%%%%%%%%%%%%%%% AEROSOL LAYER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ASLFile = 'NULL';                % ASLFile: Wavelength, SSA, a1, a+, a-, a4, b1, b2
        TauMol = round(0.65*TauRay,7);   %>= 0
        TauAbs = 0;                      %>= 0
        TauSol = 0;                      %>= 0
        SSASol = 0;                      %SSAsol = [0..1], 0 - default (from ASLFile)
        Delta = round(getAtmosphericDepolarizationFactor(wav(x)),7); %Delta = [0 , 0.9] -  depolarization factor of the atmosphere - from Bucholtz (1995)
        TKLayer = Atm_Temp;         %TKLayer = [100..6000] - Temperature of layer [Kelvin]
        fprintf(IFILE, '%-15s%-30s%-10s%-10s%-10s%-10s%-10s%-10s\n', '//AEROSOL:', 'ASLFile', 'TauMol', 'TauAbs', 'TauSol', 'SSASol', 'Delta', 'TKLayer');
        fprintf(IFILE, '%-15s%-30s%-9.7g %-9.7g %-9.7g %-9.7g %-9.7g %-9.7g \n\n', '  AEROSOL:', ASLFile, TauMol, TauAbs, TauSol, SSASol, Delta, TKLayer); 

        %%%%%%%%%%%%%%%%%%%%%%%%%% AEROSOL LAYER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ASLFile = 'oceanic.sol';    % ASLFile: Wavelength, SSA, a1, a+, a-, a4, b1, b2
        TauMol = round(0.34*TauRay,7);       %>= 0
        TauAbs = 0;                 %>= 0
        TauSol = round(0.34*(0.1*(wav(x)/440)^-1.4),7); %aot_869_modis*(WAVELEN/869)^-angstrom_modis;          %>= 0
        SSASol = SSA;               %SSAsol = [0..1], 0 - default (from ASLFile)
        Delta = 0.0279;             %Delta = [0 , 0.9] -  depolarization factor of the atmosphere 
        TKLayer = Atm_Temp;         %TKLayer = [100..6000] - Temperature of layer [Kelvin]
        fprintf(IFILE, '%-15s%-30s%-10s%-10s%-10s%-10s%-10s%-10s\n', '//AEROSOL:', 'ASLFile', 'TauMol', 'TauAbs', 'TauSol', 'SSASol', 'Delta', 'TKLayer');
        fprintf(IFILE, '%-15s%-30s%-9.7g %-9.7g %-9.7g %-9.7g %-9.7g %-9.7g \n\n', '  AEROSOL:', ASLFile, TauMol, TauAbs, TauSol, SSASol, Delta, TKLayer); 

         %%%%%%%%%%%%%%%%%%%%%%%%%%%% RECIEVER LAYER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RCRFile = 'NULL';    % RCRFile: *.amu(Tables Mu[sizeMu]):*.azm (Fi[sizeFi]), NULL(no file)
        ELEM = 11;           % ELEM = [11, 12, .. 34, 44] - (if bad, default ELEM = 11)
        Mu0 = Sun_Pol;            % acos(Mu0) = [0.0 .. 88.854008] - Sun polar angle, deg 
        Mu1 = Rec_Pol;            % acos(Mu1) = [0.0 .. 88.854008] - receiver polar angle, deg
        Fi = Rec_Az;              % Fi = [0.0 .. 360.0] - receiver azimuth, deg 
        Par4 = 0;            % Par[4] - reserved
        Par5 = 0;            % Par[5] - reserved
        fprintf(IFILE, '%-15s%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n', '//RECEIVER:', 'RCRFile', 'ELEM', 'Mu0', 'Mu1', 'Fi', 'NotUsed', 'NotUsed');
        fprintf(IFILE, '%-15s%-10s%-9.7g %-9.7g %-9.7g %-9.7g %-9.7g %-9.7g \n\n', '  RECEIVER:', RCRFile, ELEM, Mu0, Mu1, Fi, Par4, Par5);

         %%%%%%%%%%%%%%%%%%%%%%%%%% AEROSOL LAYER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ASLFile = 'oceanic.sol';         % ASLFile: Wavelength, SSA, a1, a+, a-, a4, b1, b2
        TauMol = round(0.01*TauRay,7);   %>= 0
        TauAbs = 0;                 %>= 0
        TauSol = 0;                 %aot_869_modis*(WAVELEN/869)^-angstrom_modis;          %>= 0
        SSASol = SSA;               %SSAsol = [0..1], 0 - default (from ASLFile)
        Delta = 0.0279;             %Delta = [0 , 0.9] -  depolarization factor of the atmosphere
        TKLayer = Atm_Temp;         %TKLayer = [100..6000] - Temperature of layer [Kelvin]
        fprintf(IFILE, '%-15s%-30s%-10s%-10s%-10s%-10s%-10s%-10s\n', '//AEROSOL:', 'ASLFile', 'TauMol', 'TauAbs', 'TauSol', 'SSASol', 'Delta', 'TKLayer');
        fprintf(IFILE, '%-15s%-30s%-9.7g %-9.7g %-9.7g %-9.7g %-9.7g %-9.7g \n\n', '  AEROSOL:', ASLFile, TauMol, TauAbs, TauSol, SSASol, Delta, TKLayer); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%% SURFACE LAYER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SRFtype = waveName;      % SRFtype:[NULL,ISTROPIC,DIRECTED]
                                        % NULL - smooth (all pars ignored and set to 0!)
                                        % ISOTROPIC - isotropic waves
                                        % DIRECTED - directed waves in SA part

        Wind = WS(zeta);                   % Wind, m/s = 0 - no Windy waves [2.5 .. 15] Windy waves
        WAzimuth = 0;               % WAzimuth, deg = [0 .. 360] -  Windy waves azimuth (SA only)
        TSwell = 0;                 % TSwell, s  = 0 - No Swell, [0.5 .. 15] -  Swell waves period (SA only)
        HSwell = 0;                 % HSwell, m  = [0.01 .. 5] Swell waves mean height (SA only)
        SwellAzimuth = 0;           % SwellAzimuth, deg = [0..360] - Swell waves azimuth (SA only)
        TKLayer = 291.11;              % TKLayer = [100..6000] -  Temperature of layer (on the Kelvin scale)
        fprintf(IFILE, '%-15s%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n', '//SURFACE:', 'SRFtype', 'Wind', 'WAzimuth', 'TSwell', 'HSwell', 'SwellAzi', 'TKLayer');
        fprintf(IFILE, '%-15s%-10s%-9.7g %-9.7g %-9.7g %-9.7g %-9.7g %-9.7g \n\n', '  SURFACE:', SRFtype, Wind, WAzimuth, TSwell, HSwell, SwellAzimuth, TKLayer);

        %%%%%%%%%%%%%%%%%%%%%%%% HYDROSOL LAYER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        HSLFile = HSL_str; %sim_string;  %HSLFile:  Wavelength, SSA, a1, a+, a-, a4, b1, b2
        TauMol = round(TauMol_hsol(x),7);         % TauMol >= 0
        TauAbs = round(TauAbs_hsol(x),7);         % TauAbs >= 0
        TauSol = round(TauSol_hsol(x),7);         % TauSol >= 0
        SSASol = round(SSA_hsol(x),7);         % SSAsol = [0..1], 0 - default (from HSLFile)
        Delta = round(Delta_hsol,7);          % Delta = [0..0.9] - depolarization factor for MolPart!
        Par5 = 0;           % Par[5] - reserved
        fprintf(IFILE, '%-15s%-20s%-10s%-10s%-10s%-10s%-10s%-10s\n', '//HYDROSOL:', 'HSLFile', 'TauMol', 'TauAbs', 'TauSol', 'SSASol', 'Delta', 'NotUsed');
        fprintf(IFILE, '%-15s%-20s%-9.7g %-9.7g %-9.7g %-9.7g %-9.7g %-9.7g \n\n',...
                       '  HYDROSOL:', HSLFile, TauMol, TauAbs, TauSol, SSASol, Delta, Par5); 

        %%%%%%%%%%%%%%%%%%%%%%%%%% BOTTOM LAYER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        BTMFile = 'NULL';       % BTMFile:  R[m][mu2][mu02] - MAX_NODE2
        Albedo = 0;             % Albedo = [0 .. 1]
        DiffSpec = 0;           % DiffSpec   = [0 .. 1] (Parameter for specular matrix generation)
        Spec = 0;               % Spec   = [0 .. 1] (SA only)
        VarSpec = 0;            % VarSpec = [0 .. 1] (SA only)
        Par4 = 0;               % Par[4] - reserved
        Par5 = 0;               % Par[5] - reserved
        fprintf(IFILE, '%-15s%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n', '//BOTTOM:', 'BTMFile', 'Albedo', 'DiffSpec', 'Spec', 'VarSpec', 'NotUsed', 'NotUsed');
        fprintf(IFILE, '%-15s%-10s%-9.7g %-9.7g %-9.7g %-9.7g %-9.7g %-9.7g \n\n', '  BOTTOM:', BTMFile, Albedo, DiffSpec, Spec, VarSpec, Par4, Par5);


        %%%%%%%%%% NO MORE LAYER DEFINITIONS AFTER THIS LINE %%%%%%%%%%%%%%%%%%%%%%
        %Write layer footer
        fprintf(IFILE, 'END_LAYERS:\n');

        %Write outfile section
        fprintf(IFILE, 'OUTFILE:\t%s\t%s\t%s\t%s\t%f\n', strcat('outfiles/',filename, '.out'), sun_pol_file, rec_pol_file, rec_azm_file, InclDirectSolarBeam);

        %Close the input file
        fclose(IFILE);

    end

    %% Text File of Specified Inputs

    filepath = strcat(sim_fol_name,filesep,inp_fol_name, filesep, 'Test_Inputs.txt');
    TFile = fopen(filepath, 'wt');

    fprintf(TFile, '%-26s%-2g%7s\n','Chlorophyl Concentration:',CHL,'mg/m^3');
    fprintf(TFile, 'Wavelength Range: %-3g%3s% 3g%3s\n',wav(1),'to',wav(m),'nm');
    fprintf(TFile, 'Number of Trials: %4g\n\n',m);

    if (windTag == true)
        fprintf(TFile, '%-24s%-9s\n','Waves:','Isotropic');
        fprintf(TFile, '%-24s%-4.2g%4s\n\n','Wind Speed:',WS(zeta),'m/s');
    else
        fprintf(TFile, '%-24s%-9s\n\n','Waves:','None');
        fprintf(TFile, '%-24s%-4.2g%4s\n\n','Wind Speed:',WS(zeta),'m/s');
    end

    fprintf(TFile, '%-24s%-4.2g%8s\n','Sun Polar Angle:',Sun_Pol,'degrees');
    fprintf(TFile, '%-24s%-4.2g%8s\n','Receiver Polar Angle:',Rec_Pol,'degrees');
    fprintf(TFile, '%-24s%-4.2g%8s\n\n','Receiver Azimuth Angle:',Rec_Az,'degrees');

    fprintf(TFile, 'SSA for ocean atmosphere: %4g\n\n', SSA);



    for x = 1:m

        RefractiveIndex = getRefractiveIndex(wav(x), 35, 19);      %RefractiveIndex

        fprintf(TFile, '%40s\n', '\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\');
        fprintf(TFile, 'Trial: %-3g%3s%4g\n',x,'of',m);
        fprintf(TFile, 'Wavelength: %-4g%3s\n', wav(x),'nm');     %Wavelength in nanometers
        fprintf(TFile, 'RefractiveIndex: %f\n\n', RefractiveIndex);

        fprintf(TFile, 'TauMol: %-9.7g\n',TauMol_hsol(x));
        fprintf(TFile, 'TauSol: %-9.7g\n',TauSol_hsol(x));
        fprintf(TFile, 'TauAbs: %-9.7g\n',TauAbs_hsol(x));
        fprintf(TFile, 'SSASol: %-9.7g\n\n',SSA_hsol(x));

        fprintf(TFile, 'Absorption Coefficient, a: %-9.7g\n',A_total(x));
        fprintf(TFile, 'Scattering Coefficient, b: %-9.7g\n',B_total(x));
        fprintf(TFile, 'Backscatter Coefficient, bb: %-9.7g\n\n',Bb_total(x));
    end

    fclose(TFile);

    %% Plotting
    
    if (plotTag == true)
        %Absorption
        figure(fignum)
        hold on
        subplot(1,3,1)
        plot(wav,A_total,wav,a_water,wav,a_phy,wav,a_CDOM,wav,a_NAP)
        grid on
        str1 = strcat('Absorption -',32,num2str(CHL),' mg/m^3 [Chl]');
        title(str1)
        xlabel('Wavelength (microns)')
        ylabel('Absorption Coeff. m-1')
        legend('Total','Water','Phytoplankton','CDOM','NAP')
        xlim([wav(1),wav(m)])

        %Scattering
        grid on
        subplot(1,3,2)
        plot(wav, B_total,wav,b_water,wav,b_phy,wav,b_NAP)
        str2 = strcat('Scattering -',32,num2str(CHL),' mg/m^3 [Chl]');
        title(str2)
        xlabel('Wavelength (microns)')
        ylabel('Scattering Coeff. m-1')
        legend('Total','Water','Phytoplankton','NAP')
        xlim([wav(1),wav(m)])


        %Backscatter
        grid on
        subplot(1,3,3)
        plot(wav, Bb_total,wav,bb_water,wav,bb_phy,wav,bb_NAP)
        str2 = strcat('Backscattering -',32,num2str(CHL),' mg/m^3 [Chl]');
        title(str2)
        xlabel('Wavelength (microns)')
        ylabel('Backscatter Coeff. m-1')
        legend('Total','Water','Phytoplankton','NAP')
        xlim([wav(1),wav(m)])

        fignum = fignum +1;
        
        %Rrs - Not sure what this is
        figure(fignum);
        plot(wav, 0.089*Bb_total./(A_total+Bb_total));
    end
end
clearvars