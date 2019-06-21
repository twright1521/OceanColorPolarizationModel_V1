%% |||||||||||||||DISCRIPTION |||||||||||||||||||||||||||||||||||||||||||||
%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% This program plots the data stored in 'Sim_Data.mat' 
% from the folder 'Simulation_X'

% This program is capable of producing up to 12 different plots
% of the available variables versus wavelength.
% The variables must be specified at the start of the program and are:

% Intensity, Q, U, V, Polarization, TauMol, TauSol,
% TauAbs, SSASol, A, B, and Bb

% A, B, and Bb are the absorption, scattering, and backscattering 
% coefficients, respectively

clearvars

%% |||||||||||||||INFO TO CHANGE AT THE START OF EVERY SIMULATION |||||||||
%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

sim_num = 5;          %Simulation Folder Number

% Variables to Plot

plot_vars = {'Intensity','Q','U','V','Polarization'};

[~,n] = size(plot_vars);

%% |||||||||||||||CONSTANTS |||||||||||||||||||||||||||||||||||||||||||||||
%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% Variable Names

var_names = {'[CHL]', 'Wind Speed', 'Sun Zenith', 'Receiver Zenith',...
            'Receiver Azimuth','SSA Ocean Atm','Wavelength', 'Intensity','Q',...
            'U','V','Polarization' 'TauMol', 'TauSol',...
            'TauAbs', 'SSASol', 'A', 'B', 'Bb'};
        
%% ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  P R O G R A M  S T A R T %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sim_fol_name = strcat('Simulation_',num2str(sim_num));

%Load simulation data
load(strcat(sim_fol_name,filesep,'Sim_Data.mat'))

[m,~] = size(Sim_Data_Cell);

% Make Figures Folder
 if (exist(strcat(sim_fol_name,filesep,'Figures'), 'dir') == 0)
        mkdir(strcat(sim_fol_name,filesep,'Figures'))
 end
    
% Create Legend and Half of Title
y = 1;
z = 1;
for x = 1:6
    
    for w = 1:m
        temp_array(w) = Sim_Data_Cell{w,x};
    end
    
    if (temp_array == temp_array(1))
        title_array(z) = temp_array(1);
        title_name{z} = var_names{x};
        z = z + 1;
    else
        legend_array(y,:) = temp_array;
        legend_name{y} = var_names{x};
        y = y + 1;
    end
    
    clearvars temp_array
end


% Title
title_half = 'Constants -';

for x = 1:z-1
    
    if (mod(x,2) && x ~= 1)
        title_half = sprintf('%s\n%s',title_half,...
            strcat(title_name{x},':',32,num2str(title_array(x))));
    else
        title_half = strcat(title_half,32,title_name{x},':',32,num2str(title_array(x)));
    end
end

% Legend

legend_cell = cell(1,m);
    
for w = 1:m
    
    for x = 1:y-1

        if (mod(x,2) && x ~= 1)
            legend_cell{w} = sprintf('%s\n%s',legend_cell{w},...
                strcat(legend_name{x},':',32,num2str(legend_array(x,w))));
        else
            legend_cell{w} = strcat(legend_cell{w},32,legend_name{x},':',32,num2str(legend_array(x,w)));
        end
    end
end

clearvars y z
%% |||||||||||||||PLOTTING ||||||||||||||||||||||||||||||||||||||||||||||||
fignum = 1;

for x = 1:n
    
    title_full = sprintf('%s vs. %s\n%s',plot_vars{x},var_names{7},title_half);
    
    idx = find(contains(var_names, plot_vars{x}))-6;
    
    figure(fignum)
    hold on 
    grid on
    
    for w = 1:m
        temp_array = Sim_Data_Cell{w,7};
        
        [~,k] = size(temp_array);
        
        for y = 1:k
            wav_array(y) = temp_array{1,y};
            var_array(y) = temp_array{idx,y};
        end
        
        plot(wav_array,var_array,'LineWidth',0.8)
        clearvars temparray
    end
    
    xlabel('Wavelength - nm')
    
    if (contains('Intensity',plot_vars{x}))
        ylabel(strcat(plot_vars{x},32,'- W/(m^2*um*Sr)'))
    else
        if (contains('Polarization',plot_vars{x}))
            ylabel(strcat('Degree of',32,plot_vars{x},32,'W/(m^2*um*Sr)'))
        else
            if (contains('Q', plot_vars{x}) || ...
                    contains('U', plot_vars{x}) || ...
                    contains('V', plot_vars{x}))
                ylabel(strcat(plot_vars{x},32,'Parameter - W/(m^2*um*Sr)'))
            else
                if (contains('A', plot_vars{x}) || ...
                        contains('B', plot_vars{x}) || ...
                        contains('Bb', plot_vars{x}))
                    ylabel(strcat(plot_vars{x},32,'Coefficient'))
                else
                    ylabel(plot_vars{x})
                end
            end
        end
    end
    
    title(title_full)
    legend(legend_cell)
    
    fignum = fignum + 1;
    
    savefig(strcat(sim_fol_name,filesep,'Figures',filesep,plot_vars{x}))
end