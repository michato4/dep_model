function [res] = E600(E,phi,theta,psi)
 res=sin(theta).^6.*E(1,1,7)+cos(theta).*((-6).*sin(psi).*sin(theta).^5.*E(1,2,6)+15.*cos(theta).*sin(psi).^2.*sin(theta).^4.*E(1,3,5)+(-20).*cos(theta).^2.*sin(psi).^3.*sin(theta).^3.*E(1,4,4)+15.*cos(theta).^3.*sin(psi).^4.*sin(theta).^2.*E(1,5,3)+(-6).*cos(theta).^4.*sin(psi).^5.*sin(theta).*E(1,6,2)+cos(theta).^5.*sin(psi).^6.*E(1,7,1)+(-6).*cos(psi).*sin(theta).^5.*E(2,1,6)+15.*cos(theta).*sin(2.*psi).*sin(theta).^4.*E(2,2,5)+(-60).*cos(psi).*cos(theta).^2.*sin(psi).^2.*sin(theta).^3.*E(2,3,4)+60.*cos(psi).*cos(theta).^3.*sin(psi).^3.*sin(theta).^2.*E(2,4,3)+(-30).*cos(psi).*cos(theta).^4.*sin(psi).^4.*sin(theta).*E(2,5,2)+6.*cos(psi).*cos(theta).^5.*sin(psi).^5.*E(2,6,1)+15.*cos(psi).^2.*cos(theta).*sin(theta).^4.*E(3,1,5)+(-60).*cos(psi).^2.*cos(theta).^2.*sin(psi).*sin(theta).^3.*E(3,2,4)+90.*cos(psi).^2.*cos(theta).^3.*sin(psi).^2.*sin(theta).^2.*E(3,3,3)+(-60).*cos(psi).^2.*cos(theta).^4.*sin(psi).^3.*sin(theta).*E(3,4,2)+15.*cos(psi).^2.*cos(theta).^5.*sin(psi).^4.*E(3,5,1)+(-20).*cos(psi).^3.*cos(theta).^2.*sin(theta).^3.*E(4,1,4)+60.*cos(psi).^3.*cos(theta).^3.*sin(psi).*sin(theta).^2.*E(4,2,3)+(-60).*cos(psi).^3.*cos(theta).^4.*sin(psi).^2.*sin(theta).*E(4,3,2)+20.*cos(psi).^3.*cos(theta).^5.*sin(psi).^3.*E(4,4,1)+15.*cos(psi).^4.*cos(theta).^3.*sin(theta).^2.*E(5,1,3)+(-30).*cos(psi).^4.*cos(theta).^4.*sin(psi).*sin(theta).*E(5,2,2)+15.*cos(psi).^4.*cos(theta).^5.*sin(psi).^2.*E(5,3,1)+(-6).*cos(psi).^5.*cos(theta).^4.*sin(theta).*E(6,1,2)+6.*cos(psi).^5.*cos(theta).^5.*sin(psi).*E(6,2,1)+cos(psi).^6.*cos(theta).^5.*E(7,1,1));
end