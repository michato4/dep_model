function [res_table] = get_table_of_potential_derivatives_looup_table(pos, voltages, potential_lookup_data)
% for a given position in space and a potential (voltage) on electrodes 
% (possibly complex), the function returns a table (3d array)
% containing the potential, electric field and its derivatives up to some order
    xi = 0;
    yi = 0;
    zi = 0;
    xi = find(potential_lookup_data.xg(1:end-1)-pos(1)<=0,1,'last');
    yi = find(potential_lookup_data.yg(1:end-1)-pos(2)<=0,1,'last');
    zi = find(potential_lookup_data.zg(1:end-1)-pos(3)<=0,1,'last');
    table111 = 1i*ones(size(squeeze(potential_lookup_data.lookup_table(1,1,1,:,:,:,1)))); % because of coder
    table112 = table111;
    table121 = table111;
    table122 = table111;
    table211 = table111;
    table212 = table111;
    table221 = table111;
    table222 = table111;
    for k=1:numel(voltages)       
        table111 = table111+squeeze(potential_lookup_data.lookup_table(yi,xi,zi,:,:,:,k))*voltages(k);   
        table112 = table112+squeeze(potential_lookup_data.lookup_table(yi,xi,zi+1,:,:,:,k))*voltages(k); 
        table121 = table121+squeeze(potential_lookup_data.lookup_table(yi,xi+1,zi,:,:,:,k))*voltages(k); 
        table122 = table122+squeeze(potential_lookup_data.lookup_table(yi,xi+1,zi+1,:,:,:,k))*voltages(k); 
        table211 = table211+squeeze(potential_lookup_data.lookup_table(yi+1,xi,zi,:,:,:,k))*voltages(k); 
        table212 = table212+squeeze(potential_lookup_data.lookup_table(yi+1,xi,zi+1,:,:,:,k))*voltages(k); 
        table221 = table221+squeeze(potential_lookup_data.lookup_table(yi+1,xi+1,zi,:,:,:,k))*voltages(k); 
        table222 = table222+squeeze(potential_lookup_data.lookup_table(yi+1,xi+1,zi+1,:,:,:,k))*voltages(k); 
    end
    table111 = table111-1i*ones(size(squeeze(potential_lookup_data.lookup_table(1,1,1,:,:,:,1)))); % because of coder
    table112 = table112-1i*ones(size(squeeze(potential_lookup_data.lookup_table(1,1,1,:,:,:,1)))); % because of coder
    table121 = table121-1i*ones(size(squeeze(potential_lookup_data.lookup_table(1,1,1,:,:,:,1)))); % because of coder
    table122 = table122-1i*ones(size(squeeze(potential_lookup_data.lookup_table(1,1,1,:,:,:,1)))); % because of coder
    table211 = table211-1i*ones(size(squeeze(potential_lookup_data.lookup_table(1,1,1,:,:,:,1)))); % because of coder
    table212 = table212-1i*ones(size(squeeze(potential_lookup_data.lookup_table(1,1,1,:,:,:,1)))); % because of coder
    table221 = table221-1i*ones(size(squeeze(potential_lookup_data.lookup_table(1,1,1,:,:,:,1)))); % because of coder
    table222 = table222-1i*ones(size(squeeze(potential_lookup_data.lookup_table(1,1,1,:,:,:,1)))); % because of coder

    Input_images =1i*ones(2,2,2,numel(table111));
    Input_images(1,1,1,:) = table111(:);
    Input_images(1,1,2,:) = table112(:);
    Input_images(1,2,1,:) = table121(:);
    Input_images(1,2,2,:) = table122(:);
    Input_images(2,1,1,:) = table211(:);
    Input_images(2,1,2,:) = table212(:);
    Input_images(2,2,1,:) = table221(:);
    Input_images(2,2,2,:) = table222(:);
    XI = 1+(pos(1)-potential_lookup_data.xg(xi))/(potential_lookup_data.xg(xi+1)-potential_lookup_data.xg(xi));
    YI = 1+(pos(2)-potential_lookup_data.yg(yi))/(potential_lookup_data.yg(yi+1)-potential_lookup_data.yg(yi));
    ZI = 1+(pos(3)-potential_lookup_data.zg(zi))/(potential_lookup_data.zg(zi+1)-potential_lookup_data.zg(zi));
    
    Output_images_re = ones(numel(table111),1);
    Output_images_im = ones(numel(table111),1);
    for l=1:numel(table111)
        Output_images_re(l) = interp3(real(Input_images(:,:,:,l)),XI,YI,ZI,'linear',NaN);
        Output_images_im(l) = interp3(imag(Input_images(:,:,:,l)),XI,YI,ZI,'linear',NaN);
    end
    
    res_table = reshape(Output_images_re+1i*Output_images_im,size(table111));
end
