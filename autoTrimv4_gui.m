

function [myresults, myfinalresults] = autoTrimv4_gui(sgid,lastDive,totalDives, ncDirName, createcmd, atype)
%     clc;
    addpath("myfunctions\"); % uncomment to run with Matlab / comment to compile

    imp_c_pitch_arr = [];
    imp_pitch_gain_arr = [];
    imp_cvbd_arr = [];
    imp_croll_climb_arr = [];
    imp_croll_dive_arr = [];
    imp_croll_rate_climb_arr = [];
    imp_croll_rate_dive_arr = [];
    surf_min_arr = [];
    dog_arr = [];

    myresults = []; % to be displayed in gui table
 
    n = 1;
    
    for i = lastDive:-1:(lastDive-totalDives+1)
        
        % return implied values from trim    
        [myvalues] = trimv4(sgid,i,ncDirName); % (670,290,fullpath)

        % store all in individual dive array
        individual_dive = [i myvalues.c_vbd myvalues.imp_c_vbd myvalues.c_pitch round(myvalues.imp_c_pitch,0) ...
                            round(myvalues.pitch_gain,1) round(myvalues.imp_pitch_gain,1) ...
                            myvalues.c_roll_dive myvalues.c_roll_dive_imp round(myvalues.c_roll_turn_dive,0) round(myvalues.imp_c_rollnturn_dive_ave,0) ...
                            myvalues.c_roll_climb myvalues.c_roll_climb_imp round(myvalues.c_roll_turn_climb,0) round(myvalues.imp_c_rollnturn_climb_ave,0) ...
                            myvalues.sm_cc round(myvalues.cc_surf_min,0) round(myvalues.imp_sm_cc,0) ...
                            myvalues.max_buoy myvalues.dog];

        myresults(n,:) = individual_dive; % store matrix with values [[dive1 params] [dive2 params]... [diveN params]]
    
        % store values into array
        imp_c_pitch_arr(n) = myvalues.imp_c_pitch;
        imp_pitch_gain_arr(n) = myvalues.imp_pitch_gain;
        imp_cvbd_arr(n) = myvalues.imp_c_vbd;
        imp_croll_climb_arr(n) = myvalues.c_roll_climb_imp;
        imp_croll_dive_arr(n) = myvalues.c_roll_dive_imp;
        imp_croll_rate_climb_arr(n) = myvalues.c_roll_turn_climb;
        imp_croll_rate_dive_arr(n) = myvalues.c_roll_turn_dive;
        surf_min_arr(n) = myvalues.cc_surf_min;
        dog_arr(n) = myvalues.dog;
        
        
        % Keep last dive centers
        if lastDive == i
%             fprintf('inside i=lastDive');
            last_c_vbd = myvalues.c_vbd;
            last_c_pitch = myvalues.c_pitch;
            last_pitch_gain = myvalues.pitch_gain;
            last_c_roll_dive = myvalues.c_roll_dive;
            last_c_roll_climb = myvalues.c_roll_climb;
            last_sm_cc = myvalues.sm_cc;
            last_dog = myvalues.dog;
            last_max_buoy = myvalues.max_buoy;
            
        end
    
        n = n+1;
    
    end
    
    % get implied averages
    imp_c_pitch_ave = round(arr_average(imp_c_pitch_arr,atype),0);
    imp_pitch_gain_ave = round(arr_average(imp_pitch_gain_arr,atype),1);
    imp_cvbd_ave = round(arr_average(imp_cvbd_arr,atype),0);
    imp_croll_dive_ave = round(arr_average(imp_croll_dive_arr,atype),0);
    imp_croll_climb_ave = round(arr_average(imp_croll_climb_arr,atype),0);
    imp_croll_turn_dive_ave = round(arr_average(imp_croll_rate_dive_arr,atype),0);
    imp_croll_turn_climb_ave = round(arr_average(imp_croll_rate_climb_arr,atype),0);
    
    croll_dive_ave = (imp_croll_dive_ave+imp_croll_turn_dive_ave)/2;
    croll_climb_ave = (imp_croll_climb_ave+imp_croll_turn_climb_ave)/2;
    
    surf_min_ave = round(arr_average(surf_min_arr,atype),0);
    imp_sm_cc_ave = fix(surf_min_ave + 20.0);
    
    dog_ave = round(arr_average(dog_arr,atype),2);
    if dog_ave < 1.2
        imp_max_buoy = last_max_buoy+10;
    else
        imp_max_buoy = last_max_buoy;
    end
    

    
    if createcmd 
        text = strcat('$MAX_BUOY,', num2str(round(imp_max_buoy,0)), '\n$SM_CC,', num2str(round(imp_sm_cc_ave,0)), '\n$C_VBD,', num2str(round(imp_cvbd_ave,0)), '\n$C_PITCH,', num2str(round(imp_c_pitch_ave,0)) , '\n$PITCH_GAIN,' , num2str(round(imp_pitch_gain_ave,1)) , '\n$C_ROLL_DIVE,' , num2str(round(imp_croll_dive_ave,0)) , '\n$C_ROLL_CLIMB,' , num2str(round(imp_croll_climb_ave,0)));
        output_file = strcat('sg',num2str(sgid), '_cmdfile');
        create_cmdfile(output_file, text);
    end



    myfinalresults(1,:) = [0 last_c_vbd round(imp_cvbd_ave,0)...
                                last_c_pitch round(imp_c_pitch_ave,0) last_pitch_gain round(imp_pitch_gain_ave,1)...
                                last_c_roll_dive round(imp_croll_dive_ave,0) round(imp_croll_turn_dive_ave,0) round(croll_dive_ave,0)...
                                last_c_roll_climb round(imp_croll_climb_ave,0) round(imp_croll_turn_climb_ave,0) round(croll_climb_ave,0)...
                                last_sm_cc round(surf_min_ave,0) round(imp_sm_cc_ave,0)...
                                last_max_buoy round(dog_ave,2)];
    
    
    
end

% Functions ---------------------------------------------------------------
function create_cmdfile(filename, text)
    if exist(filename, 'file')
        delete(filename);
%         fprintf('File already exists...deleting');
    else
%         fprintf('File does not exist...creating');
    end

    fout = fopen(filename,'w');
    fprintf(fout, text);
    fclose(fout);
%     fprintf('\n--> %s was created with final implied values\n\n', filename);
end
    

function average = arr_average(myArray,atype)
    n = length(myArray);
    accum = 0;
    j = 0;
    % Simple Average
    if atype == "Simple"    
        for k= 1:n
            if isnan(myArray(k))
                % if NaN do not average
            else    
                accum = accum + myArray(k);
                j = j+1; % only count if not nan
            end
        end  
        average = accum/(j); 
    end
    % Weighed Average
    if atype == "Weighed"    
        for k= 1:n
            if isnan(myArray(k))
                % if NaN do not average
            else    
                accum = accum + myArray(k)*(n+1-k);
                j = j + (n+1-k); % only count if not nan
            end
        end  
        average = accum/(j); 
    end
        
end
