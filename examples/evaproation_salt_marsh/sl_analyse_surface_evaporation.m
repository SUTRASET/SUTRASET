% this script provides all the constitutive relationship to calculate actual evaporation rate from matric potential, salinity and weather conditions
% the constitutive relationships useds here are:
% 1. kelvin equation -- calculate relative humidity based on matric potential
% 2. osmotic potential -- calculate relative humidity affected by concentration
% 3. saturated vapor density -- as a function of temperature;
% 4. surface resistance -- using Van Der Griend (1994)
% 5. salt resistance -- using Fujimaki (2006)
% 6.  mass of salt per unit area (kg/m2) -- based on salt concentration, according to Shen (2018)
% 7. liquid water saturation -- based on matric potential, according to Fayer (1995) soil water retention curve. 
% This script serves a best example code to check how actual evaporation is calcuated in SUTRASET
% /home/cmzhang/goliath_macondo_all/uqiameri/Freshwaterlens_paper1/Residual_k/new_parameter_3/k1
%/home/cmzhang/goliath_macondo_all/uqiameri/Freshwaterlens_paper1/new_simulation3/K1/part6
%  pumping effect

clear all
%close all
% project name
name='PART6';


% read input files
%inp=i0('desaline');
inp  = inpObj(  name,'block_reading','Yes');
inp.get_x_nod_mtx;   % so that we can get x in a matrix form in inp.x_nod_mtx
inp.get_y_nod_mtx;   % so that we can get y in a matrix form in inp.y_nod_mtx
%%inp.get_z_nod_mtx;   % inp.z_nod_mtx
inp.get_dx_cell_mtx; % form inp.dx_cell_mtx
inp.get_dy_cell_mtx; % from inp.dy_cell_mtx
%%inp.get_dz_cell_mtx; % from inp.dz_cell_mtx this would not work as there
%%is only one z cell
top_section_area_mtx_m2=inp.dx_cell_mtx.* inp.z(1) ; % z is the width of the cell, which is always 1 in this case
volume_mtx_m3=-inp.dx_cell_mtx.* inp.z(1).*inp.dy_cell_mtx;
ET=etObj('ET');
%only parsing output 10-13
nod  = readNOD( name,'outputnumber',3,'outputfrom',50);
p_idx  = strcmp(nod(1).label,'Pressure');
c_idx  = strcmp(nod(1).label,'Concentration');
s_idx  = strcmp(nod(1).label,'Saturation');
sm_idx  = strcmp(nod(1).label,'SM');
% a group of constants. we prefer using variables rather than numbers to
% refer to constants and unit conversion, because by doing this your code
% will be understood by you in the next 10 years. 
constants=ConstantObj();

surface_node_index=-inp.iqcp( inp.iqcp<0);
%% 
pressure_surface_sutra=nod(3).terms{p_idx}(surface_node_index);
saturation_surface_sutra=nod(3).terms{s_idx}(surface_node_index);
concentration_surface_sutra=nod(3).terms{c_idx}(surface_node_index);
sm_surface_sutra=nod(3).terms{sm_idx}(surface_node_index);
matric_potential_surface_sutra=pressure_surface_sutra ...
    /constants.g/constants.rhow_pure_water;
volume_surface_cell_m3 = volume_mtx_m3 (surface_node_index);



surface_area_sutra=inp.dx_cell_mtx(surface_node_index)*inp.z(1);

surface_area_difference=surface_area_sutra-top_section_area_mtx_m2(1,:);
max_surface_area_difference=max(surface_area_difference);

% Please note that the P output in SUTRA has a unit of pascal, while psim
% has a unit of meters
saturation_surface_consrela=ET.get_saturation('psim', ...
    matric_potential_surface_sutra,'nreg',1);

% the error associated with the result from the model and the result from
% matlab constitutive relationship
saturation_error= saturation_surface_consrela - saturation_surface_sutra ; 

max_saturation_error=max(saturation_error);

surface_resistance_sPm=func.rs1994(saturation_surface_sutra,inp.por(1) );

relative_humidity_soil= ...
    func.rh_matric(matric_potential_surface_sutra, ET.tmi+constants.kelvin) .* ...
    func.rh_osmotic(concentration_surface_sutra);

%relative_humidity_soil= ...
%    func.rh_matric(matric_potential_surface_sutra, ET.tmi+constants.kelvin);


%inp.chi1

%inp.chi2

%code to check the location of differnt nregs
% figure;scatter(inp.x ( inp.nreg==1),inp.y ( inp.nreg==1),'r.');hold on;scatter(inp.x ( inp.nreg==2),inp.y ( inp.nreg==2),'b.')



cs_kgPkg=adsorption_freundlich(concentration_surface_sutra, inp.chi1, inp.chi2);

sm_surface_cell_const_rela = (1- inp.por(1) ) *inp.rhos * cs_kgPkg.* volume_surface_cell_m3';


sm_error=sm_surface_cell_const_rela-sm_surface_sutra;


salt_weight_per_area_kgPm2=sm_surface_cell_const_rela./surface_area_sutra';
salt_resistance_sPm=func.salt_resistance_fujimaki_sPm(salt_weight_per_area_kgPm2);

solid_salt_thickness_surface_m= salt_weight_per_area_kgPm2/constants.density_solid_nacl_kgPm3;

%%plot result

%%plot result
%subplot(1,2,2);plot(concentration_surface_sutra)

vapour_density_deficit_surface_kgPm3= func.sat_vapor_density_kgPm3(ET.tmi+constants.kelvin)*relative_humidity_soil - ...
    func.sat_vapor_density_kgPm3(ET.tma+constants.kelvin)*ET.rh;

total_resistances_surface_sPm=ET.ravt+surface_resistance_sPm+surface_resistance_sPm +salt_resistance_sPm;
%total_resistances_surface_sPm=ET.ravt+surface_resistance_sPm+surface_resistance_sPm;% +salt_resistance_sPm;

actual_evaporation_mPs=vapour_density_deficit_surface_kgPm3/constants.rhow_pure_water./total_resistances_surface_sPm;

bcof = readBCOF(name,'outputnumber',3,'outputfrom',50);
%bcof may contain user specified flux (e.g., fresh groundwater input from upland)
bctime_cell_mask=strfind(bcof(1).ibc,'BCTIME');

% fill empty cell with 0
empty_cell_index=cellfun('isempty',bctime_cell_mask);
bctime_cell_mask(empty_cell_index)={0};


%bctime_cell_mask(empty_index)={[0]};
evaporation_cell_index=find(cell2mat(bctime_cell_mask));

evaporation_sutra_kgPs=bcof(1).qin(evaporation_cell_index);
evaporation_sutra_mPs=-evaporation_sutra_kgPs./surface_area_sutra'/constants.rhow_pure_water;

%acutal_evaporation_sutra_mPs=bcof(1).qin/
%figure;plot(actual_evaporation_mPs*constants.ms2mmday);hold on
%plot(evaporation_sutra_mPs)
%func.rh_osmotic(0.265)


%figure;plot(bcof(1).qin)


% figure;plot(evaporation_sutra_kgPs)

figure;plot(evaporation_sutra_mPs*constants.ms2mmday);hold on
plot(actual_evaporation_mPs*constants.ms2mmday,'r')

figure;subplot(5,1,1);plot(salt_resistance_sPm);ylabel('salt resis\n (s//m)')
subplot(5,1,2);plot(concentration_surface_sutra);ylabel('conccentration')
subplot(5,1,3);plot(solid_salt_thickness_surface_m);ylabel('solid salt thickness (m)')
subplot(5,1,4);plot(saturation_surface_consrela);ylabel('surface sat.')
subplot(5,1,5);plot(evaporation_sutra_mPs*constants.ms2mmday);ylabel('ev (mm//day)');hold on
plot(actual_evaporation_mPs*constants.ms2mmday,'r')
%
%func.rh_matric([-10,-4000],[20,30]+273.15)

%func.rh_osmotic([0.265,0.035])


%HO = EXP(-WMW*2.D0*CHI(CC)*CC/STM);

%      TSK   -- TEMPERATURE AT THAT NODE
%     TS   -- SOIL TEMPERATURE [CELSIUS]
%     DP   -- VAN'T HOFF DISSOCIATION FACTOR
%     STM  -- MOLECULAR WEIGHT OF NACL [KG/MOL]
%     WMW   -- MOLECULAR WEIGHT OF WATER [KG/MOL]

%ele  = readELE( name);
%bcop = readBCOP(name);
%
%% index for p c and s
%xnod_idx  = strcmp(nod(1).label,'X');
%ynod_idx  = strcmp(nod(1).label,'Y');
%znod_idx  = strcmp(nod(1).label,'Z');



%%sm_idx = strcmp(nod(1).label,'SM');
%
%
%xele_idx = strcmp(ele(1).label,'X origin');
%yele_idx = strcmp(ele(1).label,'Y origin');
%%zele_idx = strcmp(ele(1).label,'Z origin');
%
%vx_idx = strcmp(ele(1).label,'X velocity');
%vy_idx = strcmp(ele(1).label,'Y velocity');
%%vz_idx = strcmp(ele(1).label,'Z velocity');

