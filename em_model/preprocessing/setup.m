%% construction of the lookup-table for potential derivatives
new_electrode_array = 0;
new_lookup_table = 0;
if exist(fullfile('..','potential_lookup_data.mat'), 'file') == 2
    potential_lookup_data = load(fullfile('..','potential_lookup_data.mat'));
else
    new_lookup_table = 1;
end
if new_lookup_table == 1 || ~strcmp(potential_lookup_data.electrode_array,electrode_array) || ~isequal(xg,potential_lookup_data.xg) || ~isequal(yg,potential_lookup_data.yg) || ~isequal(zg,potential_lookup_data.zg)
    fprintf('Construction of lookup-table for potential derivatives:\n');
end
if new_lookup_table == 1 || ~strcmp(potential_lookup_data.electrode_array,electrode_array)
    new_electrode_array = 1;
    fprintf('\t* running simulation for determining the exact potential boundary conditions for the given electrode array\n');
    model = mphload('exact_boundary_conditions.mph');
    model.param.set('em',em);
    model.param.set('sm',sm);
    model.geom('geom1').feature('imp1').set('filename',fullfile(pwd,'geometry',[electrode_array '.dxf']));
    model.geom('geom1').feature('imp1').importData();
    model.mesh.run;
    model.study('std1').run;
    model.result().export('data1').run();
    mphsave();

    comsol_data = load('exact_boundary_conditions.txt');
    xg_whole = linspace(-1000e-6,1000e-6,25); % this should span the whole manipulation area
    yg_whole = xg_whole;
    sz = abs(xg_whole(1)-xg_whole(2));
    for k=1:4 % loop over the 4 electrodes
        [X,Y] = meshgrid(xg_whole,yg_whole);
        % interpolate data into a grid (already gridded data)
        Xs = reshape(comsol_data(:,1),[sqrt(numel(comsol_data(:,1))),sqrt(numel(comsol_data(:,1)))]).';
        Ys = reshape(comsol_data(:,2),[sqrt(numel(comsol_data(:,1))),sqrt(numel(comsol_data(:,1)))]).';
        Zs = reshape(comsol_data(:,4),[sqrt(numel(comsol_data(:,1))),sqrt(numel(comsol_data(:,1)))]).';
        Z = interp2(Xs,Ys,Zs,X,Y);

        Z = rot90(Z,k-1); % utilize the symmetry of the quadrupolar array
        % get rid of the too low contributions (and reduce thus the size of stored data)
        mask = Z<.01*max(max(Z)); % less than 5% of max value will not be considered
        X(mask) = NaN; Y(mask) = NaN; Z(mask) = NaN;
        X(mask) = []; Y(mask) = []; Z(mask) = [];
        el(k).X = X;
        el(k).Y = Y;
        el(k).Z = Z;
    end
    save('exact_boundary_conditions.mat','el','sz');
    delete exact_boundary_conditions.txt
end
clear comsol_data X Xs Y Ys Z Zs

if new_electrode_array==1 || ~isequal(xg,potential_lookup_data.xg) || ~isequal(yg,potential_lookup_data.yg) || ~isequal(zg,potential_lookup_data.zg)
    exact_boundary_conditions = load('exact_boundary_conditions.mat');
    if new_electrode_array || exist(fullfile(cd,'get_table_of_potential_derivatives_mex'), 'file') ~= 3
%     if exist(fullfile(cd,'get_table_of_potential_derivatives_mex'), 'file') ~= 3
%         aux_struct.X = coder.typeof(1,[1 inf],[0 1]);
%         aux_struct.Y = coder.typeof(1,[1 inf],[0 1]);
%         aux_struct.Z = coder.typeof(1,[1 inf],[0 1]);
%         aux = [];
%         aux.el = coder.typeof(aux_struct,[4 1],0);
%         aux.sz = coder.typeof(1,1,0);
%         ebc = coder.typeof(aux);  

        fprintf('\t* compiling the core function used in the generation of the lookup-table\n');
        % prepare the functions inputs so that the MATLAB Coder can automatically
        % determine size of the variables
        pos = zeros(3,1);
        v = rand(4,1)+1i*rand(4,1);
        addpath('field_derivatives');
%         codegen get_table_of_potential_derivatives.m -args {pos,v,ebc}
        codegen get_table_of_potential_derivatives.m -args {pos,v,exact_boundary_conditions}
        rmpath('field_derivatives');
    end

    fprintf('\t* generating the lookup-table (this may take a long time)\n');
    % define a grid for the lookup table
    % (this determines the size of region, in which we will be able to
    % compute the dielectrophoretic force and torque)
%     exact_boundary_conditions = load('exact_boundary_conditions.mat');
    [X,Y,Z] = meshgrid(xg,yg,zg);
    num_of_els = 1; % 1 because only because of symmetry
    v = eye(4);
    lookup_table = cell(size(X));
    for k=1:numel(lookup_table)
        fprintf(['\t\t' num2str(k) '/' num2str(numel(lookup_table))]);
        lookup_table{k} = cell(num_of_els,1);
        for l=1:num_of_els
            lookup_table{k}{l} = get_table_of_potential_derivatives_mex([X(k);Y(k);Z(k)],v(num_of_els,:).',exact_boundary_conditions);
        end
    end
    fprintf('\t* finalizing the lookup-table\n');
    % generate the rest based on the mentioned symmetry
    flippedy = flip(lookup_table,1);
    switchxy = flip(permute(lookup_table,[2 1 3]),1);
    switchxy2 = flip(permute(lookup_table,[2 1 3]),2);
    for k=1:numel(lookup_table)
        fprintf(['\t\t' num2str(k) '/' num2str(numel(lookup_table))]);
        % 2nd electrode
        aux = ones(size(lookup_table{k}{1}));
        aux(:,2:2:end,:) = -aux(:,2:2:end,:);
        lookup_table{k}{2} = permute(switchxy{k}{1},[2 1 3]).*aux;
        % 3rd electrode (mirroing horizontal axis)
        aux = ones(size(lookup_table{k}{1}));
        aux(:,2:2:end,:) = -aux(:,2:2:end,:);
        lookup_table{k}{3} = flippedy{k}{1}.*aux;
    %     % 4th electrode
        aux = ones(size(lookup_table{k}{1}));
        aux(2:2:end,:,:) = -aux(2:2:end,:,:);
        lookup_table{k}{4} = permute(switchxy2{k}{1},[2 1 3]).*aux;
    end
    % convert the cell array format of lookuptable to multidimensional array
    new_tab = NaN*ones(size(lookup_table,1),size(lookup_table,2),size(lookup_table,3),size(lookup_table{1}{1},1),size(lookup_table{1}{1},2),size(lookup_table{1}{1},3),4);
    for k=1:size(lookup_table,1)
    for l=1:size(lookup_table,2)
    for m=1:size(lookup_table,3)
    for n=1:size(lookup_table{1}{1},1)
    for o=1:size(lookup_table{1}{1},2)
    for p=1:size(lookup_table{1}{1},2)
    for q=1:4 % num of electrodes
        new_tab(k,l,m,:,:,:,q) = lookup_table{k,l,m}{q};
    end
    end
    end
    end
    end
    end
    end
    lookup_table = new_tab;
    save(fullfile('..','potential_lookup_data.mat'),'lookup_table','xg','yg','zg','electrode_array');
    clear lookup_table new_tab
end
%% generation of the basis
new_basis = 0;
if exist(fullfile('..','multipoles_basis.mat'), 'file') == 2
    multipoles_basis = load(fullfile('..','multipoles_basis.mat'));
else
    new_basis = 1;
end  
if new_basis==1 || f~=multipoles_basis.f || em~=multipoles_basis.em || sm~=multipoles_basis.sm || ep~=multipoles_basis.ep || sp~=multipoles_basis.sp || ~strcmp(object,multipoles_basis.object)
    fprintf('Generation of the basis:\n');
    fprintf('\t* opening the model file\n');
    model = mphload('extraction_of_multipoles.mph');
    model.param.set('em',em);
    model.param.set('sm',sm);
    model.param.set('ep',ep);
    model.param.set('sp',sp);
    model.param.set('f',f);
    fprintf('\t* loading geometry of the object\n');
    model.geom('geom1').feature('imp1').set('filename',fullfile(pwd,'geometry',[object '.stp']));
    model.geom('geom1').feature('imp1').importData();
    model.geom.run();
    fprintf('\t* meshing the geometry\n');
    model.mesh.run();
    fprintf('\t* running the simulations\n');
    params_names = {'q10','q11','q20','q21','q22','q30','q31','q32','q33','q40','q41','q42','q43','q44','q50','q51','q52','q53','q54','q55'};
    real_extracted_basis_re = [];
    imag_extracted_basis_re = [];
    real_extracted_basis_im = [];
    imag_extracted_basis_im = [];
    for l=1:numel(params_names)
        fprintf(['\t\t' num2str(l) '/' num2str(numel(params_names))]);
        for m=1:2 % real and imaginary part
            if (l==1 || l==3 || l==6 || l==10 || l==15) && m==2, continue; end
            % set params
            for k=1:numel(params_names)
                if k==l
                    if m==1
                        val = '1'; 
                    else
                        val = 'i';
                    end
                else
                    val = '0'; 
                end
                model.param.set(params_names{k},val);
            end
            model.study('std1').run;
            for k=[5:8 14:21]
                model.result.table(['tbl' num2str(k)]).clearTableData();
            end
            for k=[1 11 3:10 23:33 13:22 34:44]
                if sum(k==[1 11 4 7 23 28 13 14 16 19 34 39])>0
                    model.result.numerical(['int' num2str(k)]).setResult();
                else
                    model.result.numerical(['int' num2str(k)]).appendResult();
                end
            end
            monopole_data_re = mphtable(model,'tbl5');
            dipole_data_re = mphtable(model,'tbl6');
            quadrupole_data_re = mphtable(model,'tbl7');
            octupole_data_re = mphtable(model,'tbl8');
            hexadecapole_data_re = mphtable(model,'tbl14');
            t32pole_data_re = mphtable(model,'tbl15');
            monopole_data_im = mphtable(model,'tbl16');
            dipole_data_im = mphtable(model,'tbl17');
            quadrupole_data_im = mphtable(model,'tbl18');
            octupole_data_im = mphtable(model,'tbl19');
            hexadecapole_data_im = mphtable(model,'tbl20');
            t32pole_data_im = mphtable(model,'tbl21');
            if m==1
                real_extracted_basis_re(l,:) = [dipole_data_re.data, quadrupole_data_re.data, octupole_data_re.data, hexadecapole_data_re.data, t32pole_data_re.data];
                imag_extracted_basis_re(l,:) = [dipole_data_im.data, quadrupole_data_im.data, octupole_data_im.data, hexadecapole_data_im.data, t32pole_data_im.data];
            else
                real_extracted_basis_im(l,:) = [dipole_data_re.data, quadrupole_data_re.data, octupole_data_re.data, hexadecapole_data_re.data, t32pole_data_re.data];
                imag_extracted_basis_im(l,:) = [dipole_data_im.data, quadrupole_data_im.data, octupole_data_im.data, hexadecapole_data_im.data, t32pole_data_im.data];
            end
        end
        save(fullfile('..','multipoles_basis.mat'),'real_extracted_basis_re','imag_extracted_basis_re','real_extracted_basis_im','imag_extracted_basis_im','em','sm','ep','sp','f','object');
    end
    mphsave();
    clear model
end
%% compile the core components of the model
addpath('..')
if exist('get_ft_mex', 'file') ~= 3
    % the compilation is always done for a specific PC and for a specific lookup tables and basis collections
    fprintf('Compilation of the model:\n');
    fprintf('\t* preparing for the compilation of the model\n');
    addpath('rotate_field','field_derivatives');

    % prepare the functions inputs so that the MATLAB Coder can automatically
    % determine size of the variables
    pos = [0; 0; 30e-6];
    orient = [0; 0; 0];
    ampl = 45;
    v = [ampl*exp(1i*0); ampl*exp(1i*pi/2); ampl*exp(1i*pi); ampl*exp(1i*3/2*pi)];
    
    aux = [];
    aux.electrode_array = coder.typeof('s',[1 inf],[0 1]);
    aux.lookup_table = coder.typeof(1,[inf inf inf 7 7 7 4],[1 1 1 0 0 0 0]);
    aux.xg = coder.typeof(1,[1 inf],[0 1]);
    aux.yg = coder.typeof(1,[1 inf],[0 1]);
    aux.zg = coder.typeof(1,[1 inf],[0 1]);
    pld = coder.typeof(aux);
    
    aux = [];
    aux.em = coder.typeof(1,1,0);
    aux.ep = coder.typeof(1,1,0);
    aux.f = coder.typeof(1,1,1);
    aux.imag_extracted_basis_im = coder.typeof(1+1i,[20 20],[0 0]);
    aux.imag_extracted_basis_re = coder.typeof(1+1i,[20 20],[0 0]);
    aux.object = coder.typeof('s',[1 inf],[0 1]);
    aux.real_extracted_basis_im = coder.typeof(1+1i,[20 20],[0 0]);
    aux.real_extracted_basis_re = coder.typeof(1+1i,[20 20],[0 0]);
    aux.sm = coder.typeof(1,1,0);
    aux.sp = coder.typeof(1,1,0);
    mb = coder.typeof(aux);

    fprintf('\t* compiling the model\n');
    codegen get_ft.m -args {pos, orient, v, mb, pld}
    movefile('get_ft_mex.mexw64',fullfile('..','get_ft_mex.mexw64'));
    rmpath('rotate_field','field_derivatives');
end
rmpath('..')