%% |||||||||||||||DISCRIPTION |||||||||||||||||||||||||||||||||||||||||||||
%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% This program extracts data for the specified inputs from the 
% RayXP output files produced from the generated input files
% from WriteRayXPinputs.m

% Output files must be generated in RayXP and put into the folder 
% 'Output Files' within their respective RayXP folder 
% of the format 'RayXP_CHL_X_WS_Y_V_Z'

clearvars

%% |||||||||||||||INFO TO CHANGE AT THE START OF EVERY SIMULATION |||||||||
%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

sim_num = 4;        %Simulation Folder Number


%% ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  P R O G R A M  S T A R T %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sim_fol_name = strcat('Simulation_',num2str(sim_num));          %Simulation Folder

B = dir(sim_fol_name);

[m,~] = size(B);

for x = 3: m              %Folder from which output files will be taken
    filename{1,x-2} = string(B(x).name);     
end


Sim_Data_Cell = cell(m-2,7);

k = 0;

for w = 1: m-2
    
    % Load Text File

    flag = strfind(filename{w},'RayXP');
    
    if (isempty(flag))
        Sim_Data_Cell(w-k,:) = [];
        k = k + 1;
        continue
    else
        FID1 = fopen(strcat(sim_fol_name,filesep,filename{w},filesep,'Test_Inputs.txt'), 'rt');

        s = textscan(FID1, '%s', 'delimiter', '\n');

        fclose(FID1);

        A = string(s{1});

        for x = length(A):-1:1

            if (A(x) == '')
                A(x,:) = [];
            end
        end

        A = cellstr(A);

        % Puts constants into the Sim Data cell

        cellnum = 1;

        for y = 1:9

            if (y == 2 || y == 4)
                continue
            end

            test_str = A{y};
            idx1 = strfind(test_str, ':');

            for x = 1:idx1
                test_str(1) = [];
            end

            if (y == 1)
                idx2 = strfind(test_str, 'mg/m^3');
            else
                if (y == 5)
                    idx2 = strfind(test_str, 'm/s');
                else 
                    if (y ==6 || y == 7 || y ==8)
                        idx2 = strfind(test_str, 'degrees');
                    end
                end
            end

            if ( y == 3)
                N = str2double(test_str);
            else
                if ( y == 9)
                    Sim_Data_Cell{w-k,cellnum} = str2double(test_str);
                    cellnum = cellnum + 1;
                else 
                    for x = length(test_str):-1:idx2
                        test_str(x) = [];
                    end

                    Sim_Data_Cell{w-k,cellnum} = str2double(test_str);
                    cellnum = cellnum + 1;
                end
            end
        end

        D = cell(13,N);

        % Putting TauMol, TauSol, etc into D cell

        for y = 1:N

            for z = 1:9

                if (z == 2)
                    continue
                end

                test_str = A{11*y + z};
                idx1 = strfind(test_str, ':');

                for x = 1:idx1
                    test_str(1) = [];
                end

                if (z == 1)
                    idx2 = strfind(test_str, 'nm');

                    for x = length(test_str):-1:idx2
                        test_str(x) = [];
                    end

                    D{1,y} = str2double(test_str);
                else

                    D{z+4,y} = str2double(test_str);
                end
            end
        end

        % Getting I,Q,U,V,P from output file and putting into D cell

        Sun_Zen = Sim_Data_Cell{w-k,3};

        Rec_Azm = Sim_Data_Cell{w-k,5};

        Rec_Zen = Sim_Data_Cell{w-k,4};

        for y = 1:N

            outfile = strcat(sim_fol_name,filesep,filename{w},filesep,'Output Files',filesep,sprintf('Test_Script_%u.out',D{1,y}));

            [Results] =  SingleValueRayXP(outfile,Sun_Zen,Rec_Azm,Rec_Zen);

            D{2,y} = Results(3);
            D{3,y} = Results(4);
            D{4,y} = Results(5);
            D{5,y} = Results(6);
            D{6,y} = Results(7);

        end

        % Puts D cell into Sim Data cell
        Sim_Data_Cell{w-k,7} = D;

        clearvars D
    end
end

Sim_Data_Cell = sortrows(Sim_Data_Cell,[1 2 3 4 5 6]);

save(strcat(sim_fol_name,filesep,'Sim_Data.mat'),'Sim_Data_Cell')

clearvars
