clear

fil=readFIL();
%basename=fil.filename(2:end-5);


ET=etObj('ET');
inp  = inpObj(fil.basename);
inp.get_x_nod_mtx;
inp.get_y_nod_mtx;
inp.get_dx_cell_mtx;
inp.get_dy_cell_mtx;
nod  = readNOD(fil.basename);
ele  = readELE(fil.basename);
bcof = readBCOF(fil.basename);
bcop = readBCOP(fil.basename);


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


save data.mat
