clear % clear all the previous simulation result so 

fil=readFIL();  % read SUTRA.FIL to find all the filenames
%basename=fil.filename(2:end-5);


ET=etObj('ET'); % read ET.inp file
inp  = inpObj(fil.basename); %read inp file
% inp  = inpObj(fil.basename,'block_reading','yes'); % this may be faster if inp file is too slow to parse
inp.get_x_nod_mtx; % this creates a inp.x_nod_mtx, which is x coordinates in matrix form
inp.get_y_nod_mtx;
inp.get_dx_cell_mtx; %this creates variable inp.x_cell_mtx, which is the length of each cell
inp.get_dy_cell_mtx; %this creates variable inp.y_cell_mtx, which is the height of each cell
nod  = readNOD(fil.basename); %read nod file
% one could also do:
%[nod,nod2]  = readNOD( name,'outputnumber',3,'outputfrom',10);
%the reading starting from the 3rd outputs and only reads up to 13th output
ele  = readELE(fil.basename);
bcof = readBCOF(fil.basename);
bcop = readBCOP(fil.basename);

% find which column in NOD file is for sm 
% sm solid salt precipitated on the soil surface (kg)
% smass solute mass in each cell (kg)
% water mass in each cell (kg)
sm_idx    = find(strcmp(nod(1).label,'SM'));
smass_idx = find(strcmp(nod(1).label,'SMASS'));
wmass_idx = find(strcmp(nod(1).label,'WMASS'));


% index for p c and s
x_nod_idx  = strcmp(nod(1).label,'X');
y_nod_idx  = strcmp(nod(1).label,'Y');
z_nod_idx  = strcmp(nod(1).label,'Z');
p_nod_idx  = strcmp(nod(1).label,'Pressure');
c_nod_idx  = strcmp(nod(1).label,'Concentration');
s_nod_idx  = strcmp(nod(1).label,'Saturation');
sm_nod_idx = strcmp(nod(1).label,'SM');


x_ele_idx = strcmp(ele(1).label,'X origin');
y_ele_idx = strcmp(ele(1).label,'Y origin');
z_ele_idx = strcmp(ele(1).label,'Z origin');

vx_ele_idx = strcmp(ele(1).label,'X velocity');
vy_ele_idx = strcmp(ele(1).label,'Y velocity');
vz_ele_idx = strcmp(ele(1).label,'Z velocity');

% check how evaporation will be affected by the change of soil saturation (suction increase)
% you may check how to use these sub commands by 
% help plot_et_range_due_to_desaturation 
% in command line
ET.plot_et_range_due_to_desaturation('nreg',1)

% save matlab environment in binary file
save data.mat
