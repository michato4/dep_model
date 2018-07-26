function [res] = E001(E,phi,theta,psi)
 res=cos(phi).*cos(theta).*E(1,1,2)+((-1).*cos(psi).*sin(phi)+cos(phi).*sin(psi).*sin(theta)).*E(1,2,1)+(sin(phi).*sin(psi)+cos(phi).*cos(psi).*sin(theta)).*E(2,1,1);
end