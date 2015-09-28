clear all
close all
% project name
name='desaline';

% read input files
%inp=i0('desaline');
inp  = inpObj(  name);
nod  = readNOD( name);
ele  = readELE( name);
bcop = readBCOP(name);
bcof = readBCOF(name);

% index for p c and s
x_idx  = strcmp(nod(1).label,'X');
y_idx  = strcmp(nod(1).label,'Y');
p_idx  = strcmp(nod(1).label,'Pressure');
c_idx  = strcmp(nod(1).label,'Concentration');
s_idx  = strcmp(nod(1).label,'Saturation');
sm_idx = strcmp(nod(1).label,'SM');
vy_idx = strcmp(ele(1).label,'Y velocity');
yele_idx = strcmp(ele(1).label,'Y origin');
