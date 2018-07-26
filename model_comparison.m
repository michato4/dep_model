%% this script includes a demo of the proposed EM model and its comparison 
%  and validation against MST numerical simulation. For details concerning 
%  usage and theoretical background of the computations,
%  please see the README.txt file.
clear all; clc;
%% set the model inputs
% position of the object (its center of gravity)
x = 0; % m
y = 0; % m
z = 140e-6; % m
% orientation of the object
psi = 10*pi/180; % rad, rotation about the x-axis
theta = 30*pi/180; % rad, rotation about the y-axis
phi = 50*pi/180; % rad, rotation about the z-axis
% voltage on electrodes
voltages = [10*exp(1i*0); 10*exp(1i*pi/2); 50*exp(1i*pi); 50*exp(1i*3/2*pi)]; % phasors defining the harmonic signals applied to the electrodes

% It is of course possible to change also the frequency of the harmonic
% voltage signal, shape of the object or its material properties, but in
% such a case it is necessary to recompute also the whole basis, which may 
% take several minutes.
f = 25e3; % Hz, frequency of the driving sinusoidal voltages
em = 24.3; % relative permittivity of medium
sm = 2.18e-4; % S/m, electric conductivity of the medium
ep = 3.2; % relative permittivity of the particle (object of interest)
sp = 5.556e-15; % S/m, electric conductivity of the particle (object of interest)
object = 'tetris_s'; % name of the file containg the object's geometry

% grid of the lookup table for electric field calculation (used then for
% interpolation)
% (This determines the range of coordinates in which we will be able to
% compute the force and torque. If you change these values, the lookup-table 
% will have to be recomputed which can take based on the number of points from 
% several hours to a few days.)
xg = -50e-6:5e-6:50e-6;
yg = xg;
zg = xg+100e-6;

% The change of the electrode array is also quite time demanding, because
% it leads also to the recomputation of the lookup-table.
electrode_array = 'quadrupolar_electrode_array'; % name of the file containg the object's geometry
%% check the inputs and in case of neccessity recalculate the basis or lookup tables for potentials; compile the model if run for the first time
% this step might take depending on the above settings from just a few minutes to units of days on a conventional PC
% check the inputs
if x<min(xg) || x>max(xg) || y<min(yg) || y>max(yg) || z<min(zg) || z>max(zg)
    error('%s\n%s','The position of the object has to lay inside the limits of the lookup table grid.',...
           ['(' num2str(min(xg)) '<=x<=' num2str(max(xg)) ', ' num2str(min(yg)) '<=y<=' num2str(max(yg)) ', ' num2str(min(zg)) '<=z<=' num2str(max(zg)) ')']);
end
% run setup.m if not done yet
cd('em_model/preprocessing');
setup;
cd('../..');
%% run the effective multipole model
% load the precomputed lookup-table and basis of multipolar moments
addpath('em_model');
potential_lookup_data = load('potential_lookup_data.mat');
multipoles_basis = load('multipoles_basis.mat');

% compute the force
t = tic;
[em_mdl_force,em_mdl_torque] = get_ft_mex([x; y; z],[psi; theta; phi],voltages,multipoles_basis,potential_lookup_data);
em_mdl_time = toc(t);

clear potential_lookup_data multipoles_basis
rmpath('em_model');
%% run the reference MST simulation
% load the model and set the parameters
model = mphload('mst_model.mph');
model.param.set('x_obj',x);
model.param.set('y_obj',y);
model.param.set('z_obj',z);
model.param.set('psi_obj',[num2str(psi) ' [rad]']);
model.param.set('theta_obj',[num2str(theta) ' [rad]']);
model.param.set('phi_obj',[num2str(phi) ' [rad]']);
model.param.set('u1',['(' num2str(voltages(1)) ') [V]']);
model.param.set('u2',['(' num2str(voltages(2)) ') [V]']);
model.param.set('u3',['(' num2str(voltages(3)) ') [V]']);
model.param.set('u4',['(' num2str(voltages(4)) ') [V]']);
model.param.set('em',em);
model.param.set('sm',sm);
model.param.set('ep',ep);
model.param.set('sp',sp);
model.param.set('f',f);
model.geom('geom1').feature('imp1').set('filename',fullfile(cd,'em_model','preprocessing','geometry',[electrode_array '.dxf']));
model.geom('geom1').feature('imp1').importData();
model.geom('geom1').feature('imp2').set('filename',fullfile(cd,'em_model','preprocessing','geometry',[object '.stp']));
model.geom('geom1').feature('imp2').importData();
model.geom.run();

% compute the force
t = tic;
model.mesh.run;
model.study('std1').run;
model.result.table('tbl1').clearTableData();
model.result.table('tbl2').clearTableData();
model.result.numerical('gev1').setResult();
model.result.numerical('gev2').appendResult();
model.result.numerical('gev3').appendResult();
model.result.numerical('int1').setResult();
model.result.numerical('int2').appendResult();
model.result.numerical('int3').appendResult();
mst_mdl_force = mphtable(model,'tbl1');
mst_mdl_force = mst_mdl_force.data(2:end); % first column represents frequency
mst_mdl_torque = mphtable(model,'tbl2');
mst_mdl_torque = mst_mdl_torque.data(2:end);
mst_mdl_time = toc(t);
%% compare the obtained results and the respective computational times
addpath('support_functions');
% compare the computational times
disp(['control-oriented model solution time: ' num2str(em_mdl_time) 's']);
disp(['MST reference solution time: ' num2str(mst_mdl_time) 's']);

% compare the resulting forces and torques
figure; clf; hold on;
title('Comparison of MST and multipolar method');
% prepare the legend
h1 = plot([0,0],[0,0],'LineStyle','--','Color',[0 0.8 0],'LineWidth',2);
h2 = plot([0,0],[0,0],'LineStyle','--','Color',[0.8 0.8 0],'LineWidth',2);
h3 = plot([0,0],[0,0],'g-','LineWidth',2);
h4 = plot([0,0],[0,0],'y-','LineWidth',2);
axis equal;
grid;

% plot the object
t = [x; y; z];
Rx = [1 0 0; 0 cos(psi) -sin(psi); 0 sin(psi) cos(psi)];
Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
R = Rz*Ry*Rx;
[vertices, faces] = stlread(fullfile('em_model','preprocessing','geometry',[object '.stl']));
vertices = (R*vertices.'+t).';
patch('Faces',faces,'Vertices',vertices,'FaceColor',[1 0 0],'EdgeColor','none','FaceAlpha',0.5);
lighting gouraud
lightangle(-45,20);
axis equal

% plot the electrodes
aux = [1 2 3 4]; % order of labels
[~,c_Poly,~,~,~] = f_LectDxf(fullfile('em_model','preprocessing','geometry',[electrode_array '.dxf']));
for k=1:size(c_Poly,1)
    text(mean(c_Poly{k}(:,1)*1e-6),mean(c_Poly{k}(:,2)*1e-6),0,[num2str(voltages(aux(k)),'%5.1d') ' V'],'HorizontalAlignment','center','FontSize',8,'Color','red');
    plot3(c_Poly{k}(:,1)*1e-6,c_Poly{k}(:,2)*1e-6,zeros(size(c_Poly{k},1),1),'r-')
end

% plot the vectors of force and torque 
normalization_coef_F = max(norm(mst_mdl_force),norm(em_mdl_force));
normalization_coef_T = max(norm(mst_mdl_torque),norm(em_mdl_torque));
sFtot = em_mdl_force/normalization_coef_F*100e-6;
sTtot = em_mdl_torque/normalization_coef_T*100e-6;
% for EM solution
plot3([0 sFtot(1)]+t(1),[0 sFtot(2)]+t(2),[0 sFtot(3)]+t(3),'g-','LineWidth',2);
plot3([0 sTtot(1)]+t(1),[0 sTtot(2)]+t(2),[0 sTtot(3)]+t(3),'y-','LineWidth',2);
sFtotcomsol = mst_mdl_force/normalization_coef_F*100e-6;
sTtotcomsol = mst_mdl_torque/normalization_coef_T*100e-6;
% for MST solution
plot3([0 sFtotcomsol(1)]+t(1),[0 sFtotcomsol(2)]+t(2),[0 sFtotcomsol(3)]+t(3),'LineStyle','--','Color',[0 0.5 0],'LineWidth',2);
plot3([0 sTtotcomsol(1)]+t(1),[0 sTtotcomsol(2)]+t(2),[0 sTtotcomsol(3)]+t(3),'LineStyle','--','Color',[0.5 0.5 0],'LineWidth',2);

% plot the numerical values for comparison
dim = [0.55 0.6 0.3 0.3];
str = {['F_{mst} = [' regexprep(num2str(mst_mdl_force,'% .2d'),' *',' ') ']'],['F_{em} = [' regexprep(num2str(em_mdl_force,' % .2d'),' *',' ') ']']};
annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','w');
dim = [0.55 0.45 0.3 0.3];
str = {['T_{mst} = [' regexprep(num2str(mst_mdl_torque,'% .2d'),' *',' ') ']'],['T_{em} = [' regexprep(num2str(em_mdl_torque,' % .2d'),' *',' ') ']']};
annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','w');

xlim(800e-6*[-1,1]);
ylim(800e-6*[-1,1]);
zlim([0,500e-6]);
view(-45,20);
xlabel('x');
ylabel('y');
legend([h1,h2,h3,h4],{'MST force','MST torque','EM force','EM torque'},'Location','northwest')
rmpath('support_functions');