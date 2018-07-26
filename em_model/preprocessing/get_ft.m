function [F,T] = get_ft(pos, orient, voltages, multipoles_basis, potential_lookup_data)
    psi = orient(1);
    theta = orient(2);
    phi = orient(3);
    e0 = 8.8542e-12; % permittivity of vacuum

%     addpath('field_derivatives_hm');
%     field_table = get_table_of_potential_derivatives_faster_hm(src_data, pos, voltages);
%     rmpath('field_derivatives_hm')
    field_table = get_table_of_potential_derivatives_looup_table(pos,voltages,potential_lookup_data);
    
    field_table_rot = rotate_field(field_table,orient);

    r = (190e-6/0.18)/2; % r = 190e-6*5/2; % size of the sphere for applying boundary conditions
    n = 1;
    sc_dip = -1/n*4*pi*e0*multipoles_basis.em*r^(2*n+1)*n/factd(2*n-1);
    n = 2;
    sc_quad = -1/n*4*pi*e0*multipoles_basis.em*r^(2*n+1)*n/factd(2*n-1);
    n = 3;
    sc_oct = -1/n*4*pi*e0*multipoles_basis.em*r^(2*n+1)*n/factd(2*n-1);
    n = 4;
    sc_hex = -1/n*4*pi*e0*multipoles_basis.em*r^(2*n+1)*n/factd(2*n-1);
    n = 5;
    sc_t32 = -1/n*4*pi*e0*multipoles_basis.em*r^(2*n+1)*n/factd(2*n-1);
    
    % dipole
    E = [field_table_rot(2,1,1); field_table_rot(1,2,1); field_table_rot(1,1,2)];
    p1 = sc_dip*E;
    % % get the pseudo cartesian moments
    p1 = stratton_to_jackson(p1);
    % get the spherical harmonic coeficients representing the source field
    % q1Ere = jackson_to_spherical(real(p1)); % weighting coefficients
    % q1Eim = jackson_to_spherical(imag(p1));
    q1Ere = conj(jackson_to_spherical(real(p1))); % weighting coefficients
    q1Eim = conj(jackson_to_spherical(imag(p1)));

    % quadrupole
    gradE = [field_table_rot(3,1,1), field_table_rot(2,2,1), field_table_rot(2,1,2);...
             field_table_rot(2,2,1), field_table_rot(1,3,1), field_table_rot(1,2,2);...
             field_table_rot(2,1,2), field_table_rot(1,2,2), field_table_rot(1,1,3)];
    p2 = sc_quad*gradE;
    % % get the pseudo cartesian moments
    p2 = stratton_to_jackson(p2);
    % get the spherical harmonic coeficients representing the source field
    % q2Ere = jackson_to_spherical(real(p2)); % weighting coefficients
    % q2Eim = jackson_to_spherical(imag(p2));
    q2Ere = conj(jackson_to_spherical(real(p2))); % weighting coefficients
    q2Eim = conj(jackson_to_spherical(imag(p2)));

    % octupole
    grad2E = 1i*ones(3,3,3);
    grad2E(:,:,1) = [field_table_rot(4,1,1), field_table_rot(3,2,1), field_table_rot(3,1,2);...
                     field_table_rot(3,2,1), field_table_rot(2,3,1), field_table_rot(2,2,2);...
                     field_table_rot(3,1,2), field_table_rot(2,2,2), field_table_rot(2,1,3)];
    grad2E(:,:,2) = [field_table_rot(3,2,1), field_table_rot(2,3,1), field_table_rot(2,2,2);...
                     field_table_rot(2,3,1), field_table_rot(1,4,1), field_table_rot(1,3,2);...
                     field_table_rot(2,2,2), field_table_rot(1,3,2), field_table_rot(1,2,3)];
    grad2E(:,:,3) = [field_table_rot(3,1,2), field_table_rot(2,2,2), field_table_rot(2,1,3);...
                     field_table_rot(2,2,2), field_table_rot(1,3,2), field_table_rot(1,2,3);...
                     field_table_rot(2,1,3), field_table_rot(1,2,3), field_table_rot(1,1,4)];
    p3 = sc_oct*grad2E;
    % % get the pseudo cartesian moments
    p3 = stratton_to_jackson(p3);
    % get the spherical harmonic coeficients representing the source field
    % q3Ere = jackson_to_spherical(real(p3)); % weighting coefficients
    % q3Eim = jackson_to_spherical(imag(p3));
    q3Ere = conj(jackson_to_spherical(real(p3))); % weighting coefficients
    q3Eim = conj(jackson_to_spherical(imag(p3)));

    % hexadecapole
    grad3E = 1i*ones(3,3,3,3);
    for k=1:numel(grad3E)
        [I,J,K,L] = ind2sub(size(grad3E),k);
        aux = [I,J,K,L];
        grad3E(k) = field_table_rot(sum(aux==1)+1,sum(aux==2)+1,sum(aux==3)+1);
    end
    p4 = sc_hex*grad3E;
    % % get the pseudo cartesian moments
    p4 = stratton_to_jackson(p4);
    % get the spherical harmonic coeficients representing the source field
    % q4Ere = jackson_to_spherical(real(p4)); % weighting coefficients
    % q4Eim = jackson_to_spherical(imag(p4));
    q4Ere = conj(jackson_to_spherical(real(p4))); % weighting coefficients
    q4Eim = conj(jackson_to_spherical(imag(p4)));

    % 32-poles
    grad4E = 1i*ones(3,3,3,3,3);
    for k=1:numel(grad4E)
        [I,J,K,L,M] = ind2sub(size(grad4E),k);
        aux = [I,J,K,L,M];
        grad4E(k) = field_table_rot(sum(aux==1)+1,sum(aux==2)+1,sum(aux==3)+1);
    end
    p5 = sc_t32*grad4E;
    % % get the pseudo cartesian moments
    p5 = stratton_to_jackson(p5);
    % get the spherical harmonic coeficients representing the source field
    % q5Ere = jackson_to_spherical(real(p5)); % weighting coefficients
    % q5Eim = jackson_to_spherical(imag(p5));
    q5Ere = conj(jackson_to_spherical(real(p5))); % weighting coefficients
    q5Eim = conj(jackson_to_spherical(imag(p5)));
    
    % get the particle multipole (for a given orientation and external source field)
    real_extracted_basis_re = multipoles_basis.real_extracted_basis_re;
    imag_extracted_basis_re = multipoles_basis.imag_extracted_basis_re;
    real_extracted_basis_im = multipoles_basis.real_extracted_basis_im;
    imag_extracted_basis_im = multipoles_basis.imag_extracted_basis_im;

    % stackup the weighting coefficients into vectors
    wcre = [q1Ere; q2Ere; q3Ere; q4Ere; q5Ere]; % weighting coefs real
    wcim = [q1Eim; q2Eim; q3Eim; q4Eim; q5Eim]; % weighting coefs real

    % weight the basis using the weighting coefficients
    qrere_r = real(wcre).'*real_extracted_basis_re; % because we do not use monopole term
    qreim_r = real(wcre).'*imag_extracted_basis_re;
    qimre_r = real(wcim).'*real_extracted_basis_re;
    qimim_r = real(wcim).'*imag_extracted_basis_re;
    qrere_i = imag(wcre).'*real_extracted_basis_im;
    qreim_i = imag(wcre).'*imag_extracted_basis_im;
    qimre_i = imag(wcim).'*real_extracted_basis_im;
    qimim_i = imag(wcim).'*imag_extracted_basis_im;

    % dipole
    p1rere_r = spherical_to_jackson(conj(qrere_r(1:2)));
    p1reim_r = spherical_to_jackson(conj(qreim_r(1:2)));
    p1imre_r = spherical_to_jackson(conj(qimre_r(1:2)));
    p1imim_r = spherical_to_jackson(conj(qimim_r(1:2)));
    p1rere_i = spherical_to_jackson(conj(qrere_i(1:2)));
    p1reim_i = spherical_to_jackson(conj(qreim_i(1:2)));
    p1imre_i = spherical_to_jackson(conj(qimre_i(1:2)));
    p1imim_i = spherical_to_jackson(conj(qimim_i(1:2)));
    p1 = ((p1rere_r+p1rere_i)+1i*(p1reim_r+p1reim_i))+1i*((p1imre_r+p1imre_i)+1i*(p1imim_r+p1imim_i));
    p1 = jackson_to_stratton(p1);
    % quadrupole
    p2rere_r = spherical_to_jackson(conj(qrere_r(3:5)));
    p2reim_r = spherical_to_jackson(conj(qreim_r(3:5)));
    p2imre_r = spherical_to_jackson(conj(qimre_r(3:5)));
    p2imim_r = spherical_to_jackson(conj(qimim_r(3:5)));
    p2rere_i = spherical_to_jackson(conj(qrere_i(3:5)));
    p2reim_i = spherical_to_jackson(conj(qreim_i(3:5)));
    p2imre_i = spherical_to_jackson(conj(qimre_i(3:5)));
    p2imim_i = spherical_to_jackson(conj(qimim_i(3:5)));
    p2 = ((p2rere_r+p2rere_i)+1i*(p2reim_r+p2reim_i))+1i*((p2imre_r+p2imre_i)+1i*(p2imim_r+p2imim_i));
    p2 = jackson_to_stratton(p2);
    % octupole
    p3rere_r = spherical_to_jackson(conj(qrere_r(6:9)));
    p3reim_r = spherical_to_jackson(conj(qreim_r(6:9)));
    p3imre_r = spherical_to_jackson(conj(qimre_r(6:9)));
    p3imim_r = spherical_to_jackson(conj(qimim_r(6:9)));
    p3rere_i = spherical_to_jackson(conj(qrere_i(6:9)));
    p3reim_i = spherical_to_jackson(conj(qreim_i(6:9)));
    p3imre_i = spherical_to_jackson(conj(qimre_i(6:9)));
    p3imim_i = spherical_to_jackson(conj(qimim_i(6:9)));
    p3 = ((p3rere_r+p3rere_i)+1i*(p3reim_r+p3reim_i))+1i*((p3imre_r+p3imre_i)+1i*(p3imim_r+p3imim_i));
    p3 = jackson_to_stratton(p3);
    % hexadecapole
    p4rere_r = spherical_to_jackson(conj(qrere_r(10:14)));
    p4reim_r = spherical_to_jackson(conj(qreim_r(10:14)));
    p4imre_r = spherical_to_jackson(conj(qimre_r(10:14)));
    p4imim_r = spherical_to_jackson(conj(qimim_r(10:14)));
    p4rere_i = spherical_to_jackson(conj(qrere_i(10:14)));
    p4reim_i = spherical_to_jackson(conj(qreim_i(10:14)));
    p4imre_i = spherical_to_jackson(conj(qimre_i(10:14)));
    p4imim_i = spherical_to_jackson(conj(qimim_i(10:14)));
    p4 = ((p4rere_r+p4rere_i)+1i*(p4reim_r+p4reim_i))+1i*((p4imre_r+p4imre_i)+1i*(p4imim_r+p4imim_i));
    p4 = jackson_to_stratton(p4);
    % 32-pole
    p5rere_r = spherical_to_jackson(conj(qrere_r(15:20)));
    p5reim_r = spherical_to_jackson(conj(qreim_r(15:20)));
    p5imre_r = spherical_to_jackson(conj(qimre_r(15:20)));
    p5imim_r = spherical_to_jackson(conj(qimim_r(15:20)));
    p5rere_i = spherical_to_jackson(conj(qrere_i(15:20)));
    p5reim_i = spherical_to_jackson(conj(qreim_i(15:20)));
    p5imre_i = spherical_to_jackson(conj(qimre_i(15:20)));
    p5imim_i = spherical_to_jackson(conj(qimim_i(15:20)));
    p5 = ((p5rere_r+p5rere_i)+1i*(p5reim_r+p5reim_i))+1i*((p5imre_r+p5imre_i)+1i*(p5imim_r+p5imim_i));
    p5 = jackson_to_stratton(p5);

    % compute DEP force and torque on an object
    aux = field_table_rot;
    cE = conj(aux);

    % compute the three components of the DEP force (Jones)
    % dipole contribution
    F1 = 1/2*real([cE(3,1,1).*p1(1)+cE(2,2,1).*p1(2)+cE(2,1,2).*p1(3),cE(2,2,1).*p1(1)+cE(1,3,1).*p1(2)+cE(1,2,2).*p1(3),cE(2,1,2).*p1(1)+cE(1,2,2).*p1(2)+cE(1,1,3).*p1(3)]);
    % quadrupole contribution
    F2 = 1/2*real(1/2*[cE(4,1,1).*p2(1,1)+2.*cE(3,2,1).*p2(1,2)+2.*cE(3,1,2).*p2(1,3)+cE(2,3,1).*p2(2,2)+2.*cE(2,2,2).*p2(2,3)+cE(2,1,3).*p2(3,3),cE(3,2,1).*p2(1,1)+2.*cE(2,3,1).*p2(1,2)+2.*cE(2,2,2).*p2(1,3)+cE(1,4,1).*p2(2,2)+2.*cE(1,3,2).*p2(2,3)+cE(1,2,3).*p2(3,3),cE(3,1,2).*p2(1,1)+2.*cE(2,2,2).*p2(1,2)+2.*cE(2,1,3).*p2(1,3)+cE(1,3,2).*p2(2,2)+2.*cE(1,2,3).*p2(2,3)+cE(1,1,4).*p2(3,3)]);
    % octupole contribution
    F3 = 1/2*real(1/6*[cE(5,1,1).*p3(1,1,1)+3.*cE(4,2,1).*p3(1,1,2)+3.*cE(4,1,2).*p3(1,1,3)+3.*cE(3,3,1).*p3(1,2,2)+6.*cE(3,2,2).*p3(1,2,3)+3.*cE(3,1,3).*p3(1,3,3)+cE(2,4,1).*p3(2,2,2)+3.*cE(2,3,2).*p3(2,2,3)+3.*cE(2,2,3).*p3(2,3,3)+cE(2,1,4).*p3(3,3,3),cE(4,2,1).*p3(1,1,1)+3.*cE(3,3,1).*p3(1,1,2)+3.*cE(3,2,2).*p3(1,1,3)+3.*cE(2,4,1).*p3(1,2,2)+6.*cE(2,3,2).*p3(1,2,3)+3.*cE(2,2,3).*p3(1,3,3)+cE(1,5,1).*p3(2,2,2)+3.*cE(1,4,2).*p3(2,2,3)+3.*cE(1,3,3).*p3(2,3,3)+cE(1,2,4).*p3(3,3,3),cE(4,1,2).*p3(1,1,1)+3.*cE(3,2,2).*p3(1,1,2)+3.*cE(3,1,3).*p3(1,1,3)+3.*cE(2,3,2).*p3(1,2,2)+6.*cE(2,2,3).*p3(1,2,3)+3.*cE(2,1,4).*p3(1,3,3)+cE(1,4,2).*p3(2,2,2)+3.*cE(1,3,3).*p3(2,2,3)+3.*cE(1,2,4).*p3(2,3,3)+cE(1,1,5).*p3(3,3,3)]);
    % hecadecapole contribution
    F4 = 1/2*real(1/24*[cE(6,1,1).*p4(1,1,1,1)+4.*cE(5,2,1).*p4(1,1,1,2)+4.*cE(5,1,2).*p4(1,1,1,3)+6.*cE(4,3,1).*p4(1,1,2,2)+12.*cE(4,2,2).*p4(1,1,2,3)+6.*cE(4,1,3).*p4(1,1,3,3)+4.*cE(3,4,1).*p4(1,2,2,2)+12.*cE(3,3,2).*p4(1,2,2,3)+12.*cE(3,2,3).*p4(1,2,3,3)+4.*cE(3,1,4).*p4(1,3,3,3)+cE(2,5,1).*p4(2,2,2,2)+4.*cE(2,4,2).*p4(2,2,2,3)+6.*cE(2,3,3).*p4(2,2,3,3)+4.*cE(2,2,4).*p4(2,3,3,3)+cE(2,1,5).*p4(3,3,3,3),cE(5,2,1).*p4(1,1,1,1)+4.*cE(4,3,1).*p4(1,1,1,2)+4.*cE(4,2,2).*p4(1,1,1,3)+6.*cE(3,4,1).*p4(1,1,2,2)+12.*cE(3,3,2).*p4(1,1,2,3)+6.*cE(3,2,3).*p4(1,1,3,3)+4.*cE(2,5,1).*p4(1,2,2,2)+12.*cE(2,4,2).*p4(1,2,2,3)+12.*cE(2,3,3).*p4(1,2,3,3)+4.*cE(2,2,4).*p4(1,3,3,3)+cE(1,6,1).*p4(2,2,2,2)+4.*cE(1,5,2).*p4(2,2,2,3)+6.*cE(1,4,3).*p4(2,2,3,3)+4.*cE(1,3,4).*p4(2,3,3,3)+cE(1,2,5).*p4(3,3,3,3),cE(5,1,2).*p4(1,1,1,1)+4.*cE(4,2,2).*p4(1,1,1,2)+4.*cE(4,1,3).*p4(1,1,1,3)+6.*cE(3,3,2).*p4(1,1,2,2)+12.*cE(3,2,3).*p4(1,1,2,3)+6.*cE(3,1,4).*p4(1,1,3,3)+4.*cE(2,4,2).*p4(1,2,2,2)+12.*cE(2,3,3).*p4(1,2,2,3)+12.*cE(2,2,4).*p4(1,2,3,3)+4.*cE(2,1,5).*p4(1,3,3,3)+cE(1,5,2).*p4(2,2,2,2)+4.*cE(1,4,3).*p4(2,2,2,3)+6.*cE(1,3,4).*p4(2,2,3,3)+4.*cE(1,2,5).*p4(2,3,3,3)+cE(1,1,6).*p4(3,3,3,3)]);
    % 32-pole contribution
    F5 = 1/2*real(1/120*[cE(7,1,1).*p5(1,1,1,1,1)+5.*cE(6,2,1).*p5(1,1,1,1,2)+5.*cE(6,1,2).*p5(1,1,1,1,3)+10.*cE(5,3,1).*p5(1,1,1,2,2)+20.*cE(5,2,2).*p5(1,1,1,2,3)+10.*cE(5,1,3).*p5(1,1,1,3,3)+10.*cE(4,4,1).*p5(1,1,2,2,2)+30.*cE(4,3,2).*p5(1,1,2,2,3)+30.*cE(4,2,3).*p5(1,1,2,3,3)+10.*cE(4,1,4).*p5(1,1,3,3,3)+5.*cE(3,5,1).*p5(1,2,2,2,2)+20.*cE(3,4,2).*p5(1,2,2,2,3)+30.*cE(3,3,3).*p5(1,2,2,3,3)+20.*cE(3,2,4).*p5(1,2,3,3,3)+5.*cE(3,1,5).*p5(1,3,3,3,3)+cE(2,6,1).*p5(2,2,2,2,2)+5.*cE(2,5,2).*p5(2,2,2,2,3)+10.*cE(2,4,3).*p5(2,2,2,3,3)+10.*cE(2,3,4).*p5(2,2,3,3,3)+5.*cE(2,2,5).*p5(2,3,3,3,3)+cE(2,1,6).*p5(3,3,3,3,3),cE(6,2,1).*p5(1,1,1,1,1)+5.*cE(5,3,1).*p5(1,1,1,1,2)+5.*cE(5,2,2).*p5(1,1,1,1,3)+10.*cE(4,4,1).*p5(1,1,1,2,2)+20.*cE(4,3,2).*p5(1,1,1,2,3)+10.*cE(4,2,3).*p5(1,1,1,3,3)+10.*cE(3,5,1).*p5(1,1,2,2,2)+30.*cE(3,4,2).*p5(1,1,2,2,3)+30.*cE(3,3,3).*p5(1,1,2,3,3)+10.*cE(3,2,4).*p5(1,1,3,3,3)+5.*cE(2,6,1).*p5(1,2,2,2,2)+20.*cE(2,5,2).*p5(1,2,2,2,3)+30.*cE(2,4,3).*p5(1,2,2,3,3)+20.*cE(2,3,4).*p5(1,2,3,3,3)+5.*cE(2,2,5).*p5(1,3,3,3,3)+cE(1,7,1).*p5(2,2,2,2,2)+5.*cE(1,6,2).*p5(2,2,2,2,3)+10.*cE(1,5,3).*p5(2,2,2,3,3)+10.*cE(1,4,4).*p5(2,2,3,3,3)+5.*cE(1,3,5).*p5(2,3,3,3,3)+cE(1,2,6).*p5(3,3,3,3,3),cE(6,1,2).*p5(1,1,1,1,1)+5.*cE(5,2,2).*p5(1,1,1,1,2)+5.*cE(5,1,3).*p5(1,1,1,1,3)+10.*cE(4,3,2).*p5(1,1,1,2,2)+20.*cE(4,2,3).*p5(1,1,1,2,3)+10.*cE(4,1,4).*p5(1,1,1,3,3)+10.*cE(3,4,2).*p5(1,1,2,2,2)+30.*cE(3,3,3).*p5(1,1,2,2,3)+30.*cE(3,2,4).*p5(1,1,2,3,3)+10.*cE(3,1,5).*p5(1,1,3,3,3)+5.*cE(2,5,2).*p5(1,2,2,2,2)+20.*cE(2,4,3).*p5(1,2,2,2,3)+30.*cE(2,3,4).*p5(1,2,2,3,3)+20.*cE(2,2,5).*p5(1,2,3,3,3)+5.*cE(2,1,6).*p5(1,3,3,3,3)+cE(1,6,2).*p5(2,2,2,2,2)+5.*cE(1,5,3).*p5(2,2,2,2,3)+10.*cE(1,4,4).*p5(2,2,2,3,3)+10.*cE(1,3,5).*p5(2,2,3,3,3)+5.*cE(1,2,6).*p5(2,3,3,3,3)+cE(1,1,7).*p5(3,3,3,3,3)]);   

    % compute the three components of the DEP torque (Jones)
    % monopole contribution (no torque)
    % dipole contribution
    T1 = 1/2*real([cE(1,1,2).*p1(2)+(-1).*cE(1,2,1).*p1(3),(-1).*cE(1,1,2).*p1(1)+cE(2,1,1).*p1(3),cE(1,2,1).*p1(1)+(-1).*cE(2,1,1).*p1(2)]);
    % quadrupole contribution
    T2 = 1/2*real([cE(2,1,2).*p2(1,2)+(-1).*cE(2,2,1).*p2(1,3)+cE(1,2,2).*p2(2,2)+(-1).*cE(1,3,1).*p2(2,3)+cE(1,1,3).*p2(2,3)+(-1).*cE(1,2,2).*p2(3,3),(-1).*cE(2,1,2).*p2(1,1)+(-1).*cE(1,2,2).*p2(1,2)+cE(3,1,1).*p2(1,3)+(-1).*cE(1,1,3).*p2(1,3)+cE(2,2,1).*p2(2,3)+cE(2,1,2).*p2(3,3),cE(2,2,1).*p2(1,1)+(-1).*cE(3,1,1).*p2(1,2)+cE(1,3,1).*p2(1,2)+cE(1,2,2).*p2(1,3)+(-1).*cE(2,2,1).*p2(2,2)+(-1).*cE(2,1,2).*p2(2,3)]);
    % octupole contribution
    T3 = 1/2*real(1/2*[cE(3,1,2).*p3(1,1,2)+(-1).*cE(3,2,1).*p3(1,1,3)+cE(2,2,2).*p3(1,2,2)+cE(2,2,2).*p3(1,2,2)+(-1).*cE(2,3,1).*p3(1,2,3)+(-1).*cE(2,3,1).*p3(1,2,3)+cE(2,1,3).*p3(1,2,3)+cE(2,1,3).*p3(1,2,3)+(-1).*cE(2,2,2).*p3(1,3,3)+(-1).*cE(2,2,2).*p3(1,3,3)+cE(1,3,2).*p3(2,2,2)+(-1).*cE(1,4,1).*p3(2,2,3)+cE(1,2,3).*p3(2,2,3)+cE(1,2,3).*p3(2,2,3)+(-1).*cE(1,3,2).*p3(2,3,3)+(-1).*cE(1,3,2).*p3(2,3,3)+cE(1,1,4).*p3(2,3,3)+(-1).*cE(1,2,3).*p3(3,3,3),(-1).*cE(3,1,2).*p3(1,1,1)+(-1).*cE(2,2,2).*p3(1,1,2)+(-1).*cE(2,2,2).*p3(1,1,2)+cE(4,1,1).*p3(1,1,3)+(-1).*cE(2,1,3).*p3(1,1,3)+(-1).*cE(2,1,3).*p3(1,1,3)+(-1).*cE(1,3,2).*p3(1,2,2)+cE(3,2,1).*p3(1,2,3)+cE(3,2,1).*p3(1,2,3)+(-1).*cE(1,2,3).*p3(1,2,3)+(-1).*cE(1,2,3).*p3(1,2,3)+cE(3,1,2).*p3(1,3,3)+cE(3,1,2).*p3(1,3,3)+(-1).*cE(1,1,4).*p3(1,3,3)+cE(2,3,1).*p3(2,2,3)+cE(2,2,2).*p3(2,3,3)+cE(2,2,2).*p3(2,3,3)+cE(2,1,3).*p3(3,3,3),cE(3,2,1).*p3(1,1,1)+(-1).*cE(4,1,1).*p3(1,1,2)+cE(2,3,1).*p3(1,1,2)+cE(2,3,1).*p3(1,1,2)+cE(2,2,2).*p3(1,1,3)+cE(2,2,2).*p3(1,1,3)+(-1).*cE(3,2,1).*p3(1,2,2)+(-1).*cE(3,2,1).*p3(1,2,2)+cE(1,4,1).*p3(1,2,2)+(-1).*cE(3,1,2).*p3(1,2,3)+(-1).*cE(3,1,2).*p3(1,2,3)+cE(1,3,2).*p3(1,2,3)+cE(1,3,2).*p3(1,2,3)+cE(1,2,3).*p3(1,3,3)+(-1).*cE(2,3,1).*p3(2,2,2)+(-1).*cE(2,2,2).*p3(2,2,3)+(-1).*cE(2,2,2).*p3(2,2,3)+(-1).*cE(2,1,3).*p3(2,3,3)]);
    % hecadecapole contribution
    T4 = 1/2*real(1/6*[cE(4,1,2).*p4(1,1,1,2)+(-1).*cE(4,2,1).*p4(1,1,1,3)+3.*cE(3,2,2).*p4(1,1,2,2)+(-3).*cE(3,3,1).*p4(1,1,2,3)+3.*cE(3,1,3).*p4(1,1,2,3)+(-3).*cE(3,2,2).*p4(1,1,3,3)+3.*cE(2,3,2).*p4(1,2,2,2)+(-3).*cE(2,4,1).*p4(1,2,2,3)+6.*cE(2,2,3).*p4(1,2,2,3)+(-6).*cE(2,3,2).*p4(1,2,3,3)+3.*cE(2,1,4).*p4(1,2,3,3)+(-3).*cE(2,2,3).*p4(1,3,3,3)+cE(1,4,2).*p4(2,2,2,2)+(-1).*cE(1,5,1).*p4(2,2,2,3)+3.*cE(1,3,3).*p4(2,2,2,3)+(-3).*cE(1,4,2).*p4(2,2,3,3)+3.*cE(1,2,4).*p4(2,2,3,3)+(-3).*cE(1,3,3).*p4(2,3,3,3)+cE(1,1,5).*p4(2,3,3,3)+(-1).*cE(1,2,4).*p4(3,3,3,3),(-1).*cE(4,1,2).*p4(1,1,1,1)+(-3).*cE(3,2,2).*p4(1,1,1,2)+cE(5,1,1).*p4(1,1,1,3)+(-3).*cE(3,1,3).*p4(1,1,1,3)+(-3).*cE(2,3,2).*p4(1,1,2,2)+3.*cE(4,2,1).*p4(1,1,2,3)+(-6).*cE(2,2,3).*p4(1,1,2,3)+3.*cE(4,1,2).*p4(1,1,3,3)+(-3).*cE(2,1,4).*p4(1,1,3,3)+(-1).*cE(1,4,2).*p4(1,2,2,2)+3.*cE(3,3,1).*p4(1,2,2,3)+(-3).*cE(1,3,3).*p4(1,2,2,3)+6.*cE(3,2,2).*p4(1,2,3,3)+(-3).*cE(1,2,4).*p4(1,2,3,3)+3.*cE(3,1,3).*p4(1,3,3,3)+(-1).*cE(1,1,5).*p4(1,3,3,3)+cE(2,4,1).*p4(2,2,2,3)+3.*cE(2,3,2).*p4(2,2,3,3)+3.*cE(2,2,3).*p4(2,3,3,3)+cE(2,1,4).*p4(3,3,3,3),cE(4,2,1).*p4(1,1,1,1)+(-1).*cE(5,1,1).*p4(1,1,1,2)+3.*cE(3,3,1).*p4(1,1,1,2)+3.*cE(3,2,2).*p4(1,1,1,3)+(-3).*cE(4,2,1).*p4(1,1,2,2)+3.*cE(2,4,1).*p4(1,1,2,2)+(-3).*cE(4,1,2).*p4(1,1,2,3)+6.*cE(2,3,2).*p4(1,1,2,3)+3.*cE(2,2,3).*p4(1,1,3,3)+(-3).*cE(3,3,1).*p4(1,2,2,2)+cE(1,5,1).*p4(1,2,2,2)+(-6).*cE(3,2,2).*p4(1,2,2,3)+3.*cE(1,4,2).*p4(1,2,2,3)+(-3).*cE(3,1,3).*p4(1,2,3,3)+3.*cE(1,3,3).*p4(1,2,3,3)+cE(1,2,4).*p4(1,3,3,3)+(-1).*cE(2,4,1).*p4(2,2,2,2)+(-3).*cE(2,3,2).*p4(2,2,2,3)+(-3).*cE(2,2,3).*p4(2,2,3,3)+(-1).*cE(2,1,4).*p4(2,3,3,3)]);
    % 32-pole contribution
    T5 = 1/2*real(1/24*[cE(5,1,2).*p5(1,1,1,1,2)+(-1).*cE(5,2,1).*p5(1,1,1,1,3)+4.*cE(4,2,2).*p5(1,1,1,2,2)+(-4).*cE(4,3,1).*p5(1,1,1,2,3)+4.*cE(4,1,3).*p5(1,1,1,2,3)+(-4).*cE(4,2,2).*p5(1,1,1,3,3)+6.*cE(3,3,2).*p5(1,1,2,2,2)+(-6).*cE(3,4,1).*p5(1,1,2,2,3)+12.*cE(3,2,3).*p5(1,1,2,2,3)+(-12).*cE(3,3,2).*p5(1,1,2,3,3)+6.*cE(3,1,4).*p5(1,1,2,3,3)+(-6).*cE(3,2,3).*p5(1,1,3,3,3)+4.*cE(2,4,2).*p5(1,2,2,2,2)+(-4).*cE(2,5,1).*p5(1,2,2,2,3)+12.*cE(2,3,3).*p5(1,2,2,2,3)+(-12).*cE(2,4,2).*p5(1,2,2,3,3)+12.*cE(2,2,4).*p5(1,2,2,3,3)+(-12).*cE(2,3,3).*p5(1,2,3,3,3)+4.*cE(2,1,5).*p5(1,2,3,3,3)+(-4).*cE(2,2,4).*p5(1,3,3,3,3)+cE(1,5,2).*p5(2,2,2,2,2)+(-1).*cE(1,6,1).*p5(2,2,2,2,3)+4.*cE(1,4,3).*p5(2,2,2,2,3)+(-4).*cE(1,5,2).*p5(2,2,2,3,3)+6.*cE(1,3,4).*p5(2,2,2,3,3)+(-6).*cE(1,4,3).*p5(2,2,3,3,3)+4.*cE(1,2,5).*p5(2,2,3,3,3)+(-4).*cE(1,3,4).*p5(2,3,3,3,3)+cE(1,1,6).*p5(2,3,3,3,3)+(-1).*cE(1,2,5).*p5(3,3,3,3,3),(-1).*cE(5,1,2).*p5(1,1,1,1,1)+(-4).*cE(4,2,2).*p5(1,1,1,1,2)+cE(6,1,1).*p5(1,1,1,1,3)+(-4).*cE(4,1,3).*p5(1,1,1,1,3)+(-6).*cE(3,3,2).*p5(1,1,1,2,2)+4.*cE(5,2,1).*p5(1,1,1,2,3)+(-12).*cE(3,2,3).*p5(1,1,1,2,3)+4.*cE(5,1,2).*p5(1,1,1,3,3)+(-6).*cE(3,1,4).*p5(1,1,1,3,3)+(-4).*cE(2,4,2).*p5(1,1,2,2,2)+6.*cE(4,3,1).*p5(1,1,2,2,3)+(-12).*cE(2,3,3).*p5(1,1,2,2,3)+12.*cE(4,2,2).*p5(1,1,2,3,3)+(-12).*cE(2,2,4).*p5(1,1,2,3,3)+6.*cE(4,1,3).*p5(1,1,3,3,3)+(-4).*cE(2,1,5).*p5(1,1,3,3,3)+(-1).*cE(1,5,2).*p5(1,2,2,2,2)+4.*cE(3,4,1).*p5(1,2,2,2,3)+(-4).*cE(1,4,3).*p5(1,2,2,2,3)+12.*cE(3,3,2).*p5(1,2,2,3,3)+(-6).*cE(1,3,4).*p5(1,2,2,3,3)+12.*cE(3,2,3).*p5(1,2,3,3,3)+(-4).*cE(1,2,5).*p5(1,2,3,3,3)+4.*cE(3,1,4).*p5(1,3,3,3,3)+(-1).*cE(1,1,6).*p5(1,3,3,3,3)+cE(2,5,1).*p5(2,2,2,2,3)+4.*cE(2,4,2).*p5(2,2,2,3,3)+6.*cE(2,3,3).*p5(2,2,3,3,3)+4.*cE(2,2,4).*p5(2,3,3,3,3)+cE(2,1,5).*p5(3,3,3,3,3),cE(5,2,1).*p5(1,1,1,1,1)+(-1).*cE(6,1,1).*p5(1,1,1,1,2)+4.*cE(4,3,1).*p5(1,1,1,1,2)+4.*cE(4,2,2).*p5(1,1,1,1,3)+(-4).*cE(5,2,1).*p5(1,1,1,2,2)+6.*cE(3,4,1).*p5(1,1,1,2,2)+(-4).*cE(5,1,2).*p5(1,1,1,2,3)+12.*cE(3,3,2).*p5(1,1,1,2,3)+6.*cE(3,2,3).*p5(1,1,1,3,3)+(-6).*cE(4,3,1).*p5(1,1,2,2,2)+4.*cE(2,5,1).*p5(1,1,2,2,2)+(-12).*cE(4,2,2).*p5(1,1,2,2,3)+12.*cE(2,4,2).*p5(1,1,2,2,3)+(-6).*cE(4,1,3).*p5(1,1,2,3,3)+12.*cE(2,3,3).*p5(1,1,2,3,3)+4.*cE(2,2,4).*p5(1,1,3,3,3)+(-4).*cE(3,4,1).*p5(1,2,2,2,2)+cE(1,6,1).*p5(1,2,2,2,2)+(-12).*cE(3,3,2).*p5(1,2,2,2,3)+4.*cE(1,5,2).*p5(1,2,2,2,3)+(-12).*cE(3,2,3).*p5(1,2,2,3,3)+6.*cE(1,4,3).*p5(1,2,2,3,3)+(-4).*cE(3,1,4).*p5(1,2,3,3,3)+4.*cE(1,3,4).*p5(1,2,3,3,3)+cE(1,2,5).*p5(1,3,3,3,3)+(-1).*cE(2,5,1).*p5(2,2,2,2,2)+(-4).*cE(2,4,2).*p5(2,2,2,2,3)+(-6).*cE(2,3,3).*p5(2,2,2,3,3)+(-4).*cE(2,2,4).*p5(2,2,3,3,3)+(-1).*cE(2,1,5).*p5(2,3,3,3,3)]);

    % transfer the result to the original reference frame
    Rx = [1 0 0; 0 cos(psi) -sin(psi); 0 sin(psi) cos(psi)];
    Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    R = Rz*Ry*Rx;
    R = R.';
    % simply rotate the two vectors
    F1 = F1*R;
    F2 = F2*R;
    F3 = F3*R;
    F4 = F4*R;
    F5 = F5*R;
    T1 = T1*R;
    T2 = T2*R;
    T3 = T3*R;
    T4 = T4*R;
    T5 = T5*R;
    
    % get the total force
    F = F1+F2+F3+F4+F5;
    % total torque
    T = T1+T2+T3+T4+T5;
end