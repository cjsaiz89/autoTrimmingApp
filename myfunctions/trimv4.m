
function [myvalues] = trimv4(sgid, irun, ncdir)
    
    
    
    id_dive = sprintf('%03d%04d',sgid,irun);
    base_file = sprintf('p%s',id_dive); % CJS File to look up pGGGDDDD in directory indicated in dd_get_dive_data()

    [loginfo,eng,results] = dd_get_dive_data_gui(id_dive,1, ncdir); % CJS GET DATA FROM NC FILE. ('GGGDDDD',1,'C:\...fullpath)
    if (isfield(results,'processing_error') || isfield(results,'skipped_profile'))
      fprintf(1,'Processing error or skipped profile\n');
      return;
    end
    % assumes loginfo, eng, and results are available from caller
    unpack_data;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CJS PITCH TRIMMING
    
    if(isfield(eng, 'pitchCtl') == 0)    
        pitch_ctl = loginfo.gc_pitch_ctl;
        pitch_ctl = [pitch_ctl(1);pitch_ctl(:);pitch_ctl(end)];
        gc_secs = loginfo.gc_st_secs;
        gc_secs = [0.0;gc_secs(:);2000000000.0];
        sg_epoch_time_s_v = sg_epoch_time;
        eng_pitch_ctl = interp1(gc_secs,pitch_ctl,sg_epoch_time_s_v,'nearest','extrap');
        pitch_control = eng_pitch_ctl;
    else
        pitch_control = eng.pitchCtl;
    end
    
    pitch = eng.pitchAng; % measued pitch
    inv_pitch_gain = 1/loginfo.PITCH_CNV; % AD counts per cm
    pitch_gain = loginfo.PITCH_GAIN;
    
    c_pitch = loginfo.C_PITCH;
    pitch_ctl_limit = 3.5; % PARAMETER max pitch control position (cm)
    ip = find(pitch_control > -pitch_ctl_limit & pitch_control < pitch_ctl_limit); % remove extreme values
    X = [ones(size(pitch_control(ip))) pitch_control(ip)];
    A = X\pitch(ip); % regress pitch against pitch_control: A(1) is intercept, A(2) is slope
    pitch_gain_imp = A(2);
    % NOTE this assumes we take the new gain completely...no fractional adjustment
    c_pitch_adjust = ( A(1)/ pitch_gain_imp )*-inv_pitch_gain; % note sign
    imp_c_pitch = c_pitch + c_pitch_adjust;

    myvalues.imp_c_pitch = imp_c_pitch;
    myvalues.imp_pitch_gain = pitch_gain_imp;
    myvalues.c_pitch = c_pitch;
    myvalues.pitch_gain = pitch_gain;
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CJS VBD TRIMMING
    
    c_vbd = loginfo.C_VBD;
    vbd_cnts_per_cc = 1.0/loginfo.VBD_CNV;
    
    if(isfield(eng, 'vbdCC') == 0)    
        vbd_ad = loginfo.gc_vbd_ad;
        vbd_ad = [vbd_ad(1);vbd_ad(:);vbd_ad(end)];
        gc_secs = loginfo.gc_st_secs;
        gc_secs = [0.0;gc_secs(:);2000000000.0];
        vbdCC_gc = (vbd_ad - loginfo.C_VBD * loginfo.VBD_CNV);
        sg_epoch_time_s_v = sg_epoch_time;
        eng_vbd_cc = interp1(gc_secs,vbdCC_gc,sg_epoch_time_s_v,'nearest','extrap');
        vbd = eng_vbd_cc;
    else
        vbd = eng.vbdCC;
    end
    
    w  = 100.*ctr1stdiffderiv(-eng.depth, sg_time); % observed w (cm/s)
    % To try avoid internal waves, etc. look a few meters vertically around the deepest part of dive
    max_depth = max(eng.depth);
    apo_i = find(eng.depth >= max_depth - 10); % PARAMETER were we turned around at the bottom 10m
    % Where was the last place during apogee where we were still descending?
    aneg_i = find(sign(w(apo_i)) < 0,1,'last');
    if isempty(aneg_i)
      fprintf(1,'Unable to find neutral apogee point?  Missing data?\n');
      return;
    else
      % CONSIDER: for better fit, keep going backward for all contiguous negative points
      apo_i = apo_i(aneg_i:end); % eliminate early bits where she was negative or bouncing
      % Where was the first place we went positive after being negative?
      apos_i = find(sign(w(apo_i)) >= 0,1,'first');
      if isempty(apos_i)
        fprintf(1,'Unable to find neutral apogee point?  Missing data?\n');
        return;
      end
      % CONSIDER: for better fit, keep going formward for all contiguous positive points (all remaining)
      apo_i = apo_i(1:apos_i);
    end
    fit = polyfit(w(apo_i),vbd(apo_i),1); % how does vbd change with w in this range?

    cc_vbd_zw = polyval(fit,[0]); 
    neutral_depth_m = mean(eng.depth(apo_i));

    c_vbd_adjust = fix(cc_vbd_zw*vbd_cnts_per_cc); % determine adjust in AD counts
    myvalues.imp_c_vbd = c_vbd + c_vbd_adjust;

    myvalues.c_vbd = c_vbd;


    % CJS ----------------------------------------------------------------------
    % TRIM ROLL BY RATE
    
    vmtime   = eng.elaps_t;
    mp = length(vmtime);
    
    hdg = eng.head;
    hdgdiff = zeros(mp,1);
    hdgdiff(2:mp) = diff(hdg);
    hdgdiff = mod(hdgdiff,360);
    in_hdg = find(hdgdiff > 180);
    hdgdiff(in_hdg) = hdgdiff(in_hdg) - 360;
    hdg_wrapped = hdg(1) + cumsum(hdgdiff);
    turn_rate = ctr1stdiffderiv(hdg_wrapped, vmtime);
    
    c_roll_dive = loginfo.C_ROLL_DIVE;
    c_roll_climb = loginfo.C_ROLL_CLIMB;
    roll_control = eng.rollCtl;
    inv_roll_cnv = 1/loginfo.ROLL_CNV; % 35.37 AD counts per degree
    
    roll_control_counts_dive = c_roll_dive + roll_control(sg_dive_i)*inv_roll_cnv;
    Xd = [ones(size(roll_control_counts_dive)) roll_control_counts_dive];
    Cd = Xd\turn_rate(sg_dive_i);
    c_roll_turn_dive = -Cd(1)/Cd(2);
    
    roll_control_counts_climb = c_roll_climb + roll_control(sg_climb_i)*inv_roll_cnv;
    Xc = [ones(size(roll_control_counts_climb)) roll_control_counts_climb];
    Cc = Xc\turn_rate(sg_climb_i);
    c_roll_turn_climb = -Cc(1)/Cc(2);
    
    myvalues.c_roll_dive = c_roll_dive;
    myvalues.c_roll_climb = c_roll_climb;
    myvalues.c_roll_turn_dive = c_roll_turn_dive;
    myvalues.c_roll_turn_climb = c_roll_turn_climb;

    % CJS ----------------------------------------------------------------------
    % TRIM ROLL

    if(isfield(eng, 'rollCtl') == 0)    
        roll_ctl = loginfo.gc_roll_ctl;
        roll_ctl = [roll_ctl(1);roll_ctl(:);roll_ctl(end)];
        gc_secs = loginfo.gc_st_secs;
        gc_secs = [0.0;gc_secs(:);2000000000.0];
        sg_epoch_time_s_v = sg_epoch_time;
        eng_roll_ctl = interp1(gc_secs,roll_ctl,sg_epoch_time_s_v,'nearest','extrap');
        roll_control = eng_roll_ctl;
    else
        roll_control = eng.rollCtl;
    end
    
    roll = eng.rollAng;

    % Roll dive
    Xd = [ones(size(roll_control(sg_dive_i))) roll_control(sg_dive_i)];
    Bd = Xd\roll(sg_dive_i);
    c_roll_dive_imp = fix(c_roll_dive - ( Bd(1)/Bd(2) )*inv_roll_cnv);
    % Roll climb
    Xc = [ones(size(roll_control(sg_climb_i))) roll_control(sg_climb_i)];
    Bc = Xc\roll(sg_climb_i);
    c_roll_climb_imp = fix(c_roll_climb - ( Bc(1)/Bc(2) )*inv_roll_cnv);

    myvalues.imp_c_rollnturn_dive_ave = (c_roll_turn_dive+c_roll_dive_imp)/2;
    myvalues.imp_c_rollnturn_climb_ave = (c_roll_turn_climb+c_roll_climb_imp)/2;
    myvalues.c_roll_dive_imp = c_roll_dive_imp;
    myvalues.c_roll_climb_imp = c_roll_climb_imp;

    % CJS ----------------------------------------------------------------------
    % TRIM SMCC

    temp_raw = results.temperature_raw;
    salin_raw = results.salinity_raw;
    density_raw = sw_dens0(salin_raw,temp_raw);
    vmdepth = eng.depth;
    zpos = -vmdepth;
    zmax = max(zpos); % shallowest point to handle yoyo dive (normally ~ -0.5 m)
    i_1m = find(zpos > (zmax - 1.0) & zpos < zmax); % PARAMETER 1meter
    i_1m = valIndex(i_1m, density_raw);
    density_1m = mean(density_raw(i_1m)); % what is density 1m lower?
    
   
    % find apogee depth within 1 meter of the deepest point
    zmin = min(zpos);
    i_apogee = find(zpos < zmin + 1); % PARAMETER density at 1 meter from max depth 
    i_apogee = valIndex(i_apogee, density_raw);
    density_apogee = mean(density_raw(i_apogee));
    
    cc_per_m3 = 1e6;
    vbd_1m = cc_per_m3*mass*(1./density_1m - 1./density_apogee); % What vbd throw do we need to go between deep and surface points?
    cc_surf_min = vbd_1m + 150; % PARAMETER 150g (or cc) to get the antenna out of the water
    
    
    myvalues.cc_surf_min = cc_surf_min;
    myvalues.imp_sm_cc = cc_surf_min + 20.0;
    myvalues.sm_cc = loginfo.SM_CC; 
    

    % CJS ----------------------------------------------------------------------
    % GPS positions - distance over ground
    gps_lat_start = deg2rad(loginfo.GPS2_lat); % dd.dddd
    gps_lon_start = deg2rad(loginfo.GPS2_lon);
    gps_lat_end = deg2rad(loginfo.GPS_lat);
    gps_lon_end = deg2rad(loginfo.GPS_lon);

    rearth = 6371.0 ; %km
    dog = acos(sin(gps_lat_start) * sin(gps_lat_end) + cos(gps_lat_start) * cos(gps_lat_end) * cos(abs(gps_lon_end-gps_lon_start)))*rearth;
    
    myvalues.dog = dog;
    myvalues.max_buoy = loginfo.MAX_BUOY;











