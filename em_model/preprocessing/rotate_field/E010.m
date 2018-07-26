function [res] = E010(E,phi,theta,psi)
 res=cos(theta).*sin(phi).*E(1,1,2)+(cos(phi).*cos(psi)+sin(phi).*sin(psi).*sin(theta)).*E(1,2,1)+((-1).*cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin(theta)).*E(2,1,1);
end