function [Results] = SingleValueRayXP(inputFileName,Sun_Zen,Rec_Azm,Rec_Zen)

% Gets a single line of data from the RayXp output for
% specified triplet of Sun Zenith, Receiver Zenith and Azimuth angles
%
% Sun_Zen - Sun Zenith angle, degrees (integer)
% Rec_Azm - Receiver Azimuth angle, degrees (integer)
% Rec_Zen - Receiver Zenith angle, degrees (integer)

IFILE = fopen(inputFileName);
while (~feof(IFILE)) 
    %Read the line
    line = fgetl(IFILE);
    
    line_test = line(2:end-49);
    
    if (strcmp(line_test,'Sun zenith angle'))
        sun_ang = str2double(strtrim(line(56:end-1)));
        if (sun_ang == Sun_Zen)
            line2 = fgetl(IFILE);
            azm_ang = str2double(strtrim(line2(56:end-1)));
            if (azm_ang == Rec_Azm)
                break
            end
        end
    end

end

flag = false;

while (flag == false)
    line3 = fgetl(IFILE);
    
    parsenumbers = textscan(line3,'%f');
    numbers = parsenumbers{1};
    if ~isempty(numbers)
        if (numbers(1) == Rec_Zen)
            flag = true;
            Results = numbers.';
        end
    end
end

fclose(IFILE);
end
