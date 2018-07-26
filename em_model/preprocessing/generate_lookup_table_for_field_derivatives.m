%% approximate exact B.C. by a set of blocks
clear all; clc;
%% load and process the exact B.C. calculated in Comsol model for later use
comsol_data = load('exact_boundary_conditions.txt');
%%
% generate a cut-out used for approximation of potential
% (the more elements, the better approximation of potential)
xg = linspace(-1000e-6,1000e-6,25); % this should span the whole manipulation area
yg = xg;
sz = abs(xg(1)-xg(2));

for k=1:4 % over the electrodes
    [X,Y] = meshgrid(xg,yg);
    % interpolate data into a grid (already gridded data)
    Xs = reshape(comsol_data(:,1),[sqrt(numel(comsol_data(:,1))),sqrt(numel(comsol_data(:,1)))]).';
    Ys = reshape(comsol_data(:,2),[sqrt(numel(comsol_data(:,1))),sqrt(numel(comsol_data(:,1)))]).';
    Zs = reshape(comsol_data(:,4),[sqrt(numel(comsol_data(:,1))),sqrt(numel(comsol_data(:,1)))]).';
    Z = interp2(Xs,Ys,Zs,X,Y);
    
    Z = rot90(Z,k-1); % take advantage of the symmetry of the electrode array
    % get rid of the too low contributions (and reduce thus the size of stored data while keeping the solution accurate enough)
    mask = Z<.01*max(max(Z)); % less than 5% of max value will not be considered
    X(mask) = NaN; Y(mask) = NaN; Z(mask) = NaN;
    % surf(X,Y,Z); pause;
    X(mask) = []; Y(mask) = []; Z(mask) = [];
    el(k).X = X;
    el(k).Y = Y;
    el(k).Z = Z;
    
    % inspect the g
    figure(1); clf;
    scatter3(X,Y,Z);
    xlim([xg(1) xg(end)]);
    ylim([yg(1) yg(end)]);    
    grid;
    pause;
end
%%
save('quadrupolar_electrode_array_bc.mat','el','sz');
%% compute the field derivatives (using analytical formula) in advance and store it in a lookup table
% define a grid for the lookup table
xg = -50e-6:5e-6:50e-6;
yg = xg;
zg = 25e-6:5e-6:125e-6;
[X,Y,Z] = meshgrid(xg,yg,zg);
%%
num_of_els = 1; % 1 because only because of symmetry
voltages = eye(4);
lookup_table = cell(size(X));
src_data = load('quadrupolar_electrode_array_bc.mat');
for k=1:numel(lookup_table)
    disp([num2str(k) '/' num2str(numel(lookup_table))]);
    lookup_table{k} = cell(num_of_els,1);
    for l=1:num_of_els
        lookup_table{k}{l} = get_table_of_potential_derivatives(src_data,[X(k);Y(k);Z(k)],voltages(num_of_els,:).');
    end
end

% generate the rest based on the mentioned symmetry
flippedx = flip(lookup_table,2);
flippedxy = flip(flippedx,1);
flippedy = flip(lookup_table,1);
for k=1:numel(lookup_table)
    disp([num2str(k) '/' num2str(numel(lookup_table))]);    
    % 2nd electrode (mirroing about vertical axis)
    aux = ones(size(lookup_table{k}{1}));
    aux(2:2:end,:,:) = -aux(2:2:end,:,:);
    lookup_table{k}{2} = flippedx{k}{1}.*aux;
    % 3rd electrode (mirroing about both horizontal and vertical axis)
    aux = ones(size(lookup_table{k}{1}));
    aux(2:2:end,:,:) = -aux(2:2:end,:,:);
    aux(:,2:2:end,:) = -aux(:,2:2:end,:);
    lookup_table{k}{3} = flippedxy{k}{1}.*aux;
    % 4th electrode (mirroing about horizonal axis)
    aux = ones(size(lookup_table{k}{1}));
    aux(:,2:2:end,:) = -aux(:,2:2:end,:);
    lookup_table{k}{4} = flippedy{k}{1}.*aux;
end
save('lookup_table_quadrupolar_electrode_array.mat','lookup_table');