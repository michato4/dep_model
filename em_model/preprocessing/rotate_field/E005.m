function [res] = E005(E,phi,theta,psi)
 res=cos(phi).^5.*cos(theta).^5.*E(1,1,6)+5.*cos(phi).^4.*cos(theta).^4.*((-1).*cos(psi).*sin(phi)+cos(phi).*sin(psi).*sin(theta)).*E(1,2,5)+10.*cos(phi).^3.*cos(psi).^2.*cos(theta).^3.*sin(phi).^2.*E(1,3,4)+(-10).*cos(phi).^4.*cos(theta).^3.*sin(phi).*sin(2.*psi).*sin(theta).*E(1,3,4)+10.*cos(phi).^5.*cos(theta).^3.*sin(psi).^2.*sin(theta).^2.*E(1,3,4)+(-10).*cos(phi).^2.*cos(psi).^3.*cos(theta).^2.*sin(phi).^3.*E(1,4,3)+30.*cos(phi).^3.*cos(psi).^2.*cos(theta).^2.*sin(phi).^2.*sin(psi).*sin(theta).*E(1,4,3)+(-30).*cos(phi).^4.*cos(psi).*cos(theta).^2.*sin(phi).*sin(psi).^2.*sin(theta).^2.*E(1,4,3)+10.*cos(phi).^5.*cos(theta).^2.*sin(psi).^3.*sin(theta).^3.*E(1,4,3)+5.*cos(phi).*cos(psi).^4.*cos(theta).*sin(phi).^4.*E(1,5,2)+30.*cos(phi).^3.*cos(psi).^2.*cos(theta).*sin(phi).^2.*sin(psi).^2.*sin(theta).^2.*E(1,5,2)+(-20).*cos(phi).^4.*cos(psi).*cos(theta).*sin(phi).*sin(psi).^3.*sin(theta).^3.*E(1,5,2)+5.*cos(phi).^5.*cos(theta).*sin(psi).^4.*sin(theta).^4.*E(1,5,2)+(-10).*cos(phi).^2.*cos(psi).^3.*sin(phi).^3.*sin(psi).*sin(2.*theta).*E(1,5,2)+(-1).*cos(psi).^5.*sin(phi).^5.*E(1,6,1)+5.*cos(phi).*cos(psi).^4.*sin(phi).^4.*sin(psi).*sin(theta).*E(1,6,1)+(-10).*cos(phi).^2.*cos(psi).^3.*sin(phi).^3.*sin(psi).^2.*sin(theta).^2.*E(1,6,1)+10.*cos(phi).^3.*cos(psi).^2.*sin(phi).^2.*sin(psi).^3.*sin(theta).^3.*E(1,6,1)+(-5).*cos(phi).^4.*cos(psi).*sin(phi).*sin(psi).^4.*sin(theta).^4.*E(1,6,1)+cos(phi).^5.*sin(psi).^5.*sin(theta).^5.*E(1,6,1)+5.*cos(phi).^4.*cos(theta).^4.*sin(phi).*sin(psi).*E(2,1,5)+5.*cos(phi).^5.*cos(psi).*cos(theta).^4.*sin(theta).*E(2,1,5)+(-10).*cos(phi).^3.*cos(theta).^3.*sin(phi).^2.*sin(2.*psi).*E(2,2,4)+(-20).*cos(phi).^4.*cos(psi).^2.*cos(theta).^3.*sin(phi).*sin(theta).*E(2,2,4)+20.*cos(phi).^4.*cos(theta).^3.*sin(phi).*sin(psi).^2.*sin(theta).*E(2,2,4)+10.*cos(phi).^5.*cos(theta).^3.*sin(2.*psi).*sin(theta).^2.*E(2,2,4)+30.*cos(phi).^2.*cos(psi).^2.*cos(theta).^2.*sin(phi).^3.*sin(psi).*E(2,3,3)+30.*cos(phi).^3.*cos(psi).^3.*cos(theta).^2.*sin(phi).^2.*sin(theta).*E(2,3,3)+(-60).*cos(phi).^3.*cos(psi).*cos(theta).^2.*sin(phi).^2.*sin(psi).^2.*sin(theta).*E(2,3,3)+(-60).*cos(phi).^4.*cos(psi).^2.*cos(theta).^2.*sin(phi).*sin(psi).*sin(theta).^2.*E(2,3,3)+30.*cos(phi).^4.*cos(theta).^2.*sin(phi).*sin(psi).^3.*sin(theta).^2.*E(2,3,3)+30.*cos(phi).^5.*cos(psi).*cos(theta).^2.*sin(psi).^2.*sin(theta).^3.*E(2,3,3)+(-20).*cos(phi).*cos(psi).^3.*cos(theta).*sin(phi).^4.*sin(psi).*E(2,4,2)+60.*cos(phi).^3.*cos(psi).^3.*cos(theta).*sin(phi).^2.*sin(psi).*sin(theta).^2.*E(2,4,2)+(-60).*cos(phi).^3.*cos(psi).*cos(theta).*sin(phi).^2.*sin(psi).^3.*sin(theta).^2.*E(2,4,2)+(-60).*cos(phi).^4.*cos(psi).^2.*cos(theta).*sin(phi).*sin(psi).^2.*sin(theta).^3.*E(2,4,2)+20.*cos(phi).^4.*cos(theta).*sin(phi).*sin(psi).^4.*sin(theta).^3.*E(2,4,2)+20.*cos(phi).^5.*cos(psi).*cos(theta).*sin(psi).^3.*sin(theta).^4.*E(2,4,2)+(-10).*cos(phi).^2.*cos(psi).^4.*sin(phi).^3.*sin(2.*theta).*E(2,4,2)+30.*cos(phi).^2.*cos(psi).^2.*sin(phi).^3.*sin(psi).^2.*sin(2.*theta).*E(2,4,2)+5.*cos(psi).^4.*sin(phi).^5.*sin(psi).*E(2,5,1)+5.*cos(phi).*cos(psi).^5.*sin(phi).^4.*sin(theta).*E(2,5,1)+(-20).*cos(phi).*cos(psi).^3.*sin(phi).^4.*sin(psi).^2.*sin(theta).*E(2,5,1)+(-20).*cos(phi).^2.*cos(psi).^4.*sin(phi).^3.*sin(psi).*sin(theta).^2.*E(2,5,1)+30.*cos(phi).^2.*cos(psi).^2.*sin(phi).^3.*sin(psi).^3.*sin(theta).^2.*E(2,5,1)+30.*cos(phi).^3.*cos(psi).^3.*sin(phi).^2.*sin(psi).^2.*sin(theta).^3.*E(2,5,1)+(-20).*cos(phi).^3.*cos(psi).*sin(phi).^2.*sin(psi).^4.*sin(theta).^3.*E(2,5,1)+(-20).*cos(phi).^4.*cos(psi).^2.*sin(phi).*sin(psi).^3.*sin(theta).^4.*E(2,5,1)+5.*cos(phi).^4.*sin(phi).*sin(psi).^5.*sin(theta).^4.*E(2,5,1)+5.*cos(phi).^5.*cos(psi).*sin(psi).^4.*sin(theta).^5.*E(2,5,1)+10.*cos(phi).^3.*cos(theta).^3.*sin(phi).^2.*sin(psi).^2.*E(3,1,4)+10.*cos(phi).^4.*cos(theta).^3.*sin(phi).*sin(2.*psi).*sin(theta).*E(3,1,4)+10.*cos(phi).^5.*cos(psi).^2.*cos(theta).^3.*sin(theta).^2.*E(3,1,4)+(-30).*cos(phi).^2.*cos(psi).*cos(theta).^2.*sin(phi).^3.*sin(psi).^2.*E(3,2,3)+(-60).*cos(phi).^3.*cos(psi).^2.*cos(theta).^2.*sin(phi).^2.*sin(psi).*sin(theta).*E(3,2,3)+30.*cos(phi).^3.*cos(theta).^2.*sin(phi).^2.*sin(psi).^3.*sin(theta).*E(3,2,3)+(-30).*cos(phi).^4.*cos(psi).^3.*cos(theta).^2.*sin(phi).*sin(theta).^2.*E(3,2,3)+60.*cos(phi).^4.*cos(psi).*cos(theta).^2.*sin(phi).*sin(psi).^2.*sin(theta).^2.*E(3,2,3)+30.*cos(phi).^5.*cos(psi).^2.*cos(theta).^2.*sin(psi).*sin(theta).^3.*E(3,2,3)+30.*cos(phi).*cos(psi).^2.*cos(theta).*sin(phi).^4.*sin(psi).^2.*E(3,3,2)+30.*cos(phi).^3.*cos(psi).^4.*cos(theta).*sin(phi).^2.*sin(theta).^2.*E(3,3,2)+(-120).*cos(phi).^3.*cos(psi).^2.*cos(theta).*sin(phi).^2.*sin(psi).^2.*sin(theta).^2.*E(3,3,2)+30.*cos(phi).^3.*cos(theta).*sin(phi).^2.*sin(psi).^4.*sin(theta).^2.*E(3,3,2)+(-60).*cos(phi).^4.*cos(psi).^3.*cos(theta).*sin(phi).*sin(psi).*sin(theta).^3.*E(3,3,2)+60.*cos(phi).^4.*cos(psi).*cos(theta).*sin(phi).*sin(psi).^3.*sin(theta).^3.*E(3,3,2)+30.*cos(phi).^5.*cos(psi).^2.*cos(theta).*sin(psi).^2.*sin(theta).^4.*E(3,3,2)+30.*cos(phi).^2.*cos(psi).^3.*sin(phi).^3.*sin(psi).*sin(2.*theta).*E(3,3,2)+(-30).*cos(phi).^2.*cos(psi).*sin(phi).^3.*sin(psi).^3.*sin(2.*theta).*E(3,3,2)+(-10).*cos(psi).^3.*sin(phi).^5.*sin(psi).^2.*E(3,4,1)+(-20).*cos(phi).*cos(psi).^4.*sin(phi).^4.*sin(psi).*sin(theta).*E(3,4,1)+30.*cos(phi).*cos(psi).^2.*sin(phi).^4.*sin(psi).^3.*sin(theta).*E(3,4,1)+(-10).*cos(phi).^2.*cos(psi).^5.*sin(phi).^3.*sin(theta).^2.*E(3,4,1)+60.*cos(phi).^2.*cos(psi).^3.*sin(phi).^3.*sin(psi).^2.*sin(theta).^2.*E(3,4,1)+(-30).*cos(phi).^2.*cos(psi).*sin(phi).^3.*sin(psi).^4.*sin(theta).^2.*E(3,4,1)+30.*cos(phi).^3.*cos(psi).^4.*sin(phi).^2.*sin(psi).*sin(theta).^3.*E(3,4,1)+(-60).*cos(phi).^3.*cos(psi).^2.*sin(phi).^2.*sin(psi).^3.*sin(theta).^3.*E(3,4,1)+10.*cos(phi).^3.*sin(phi).^2.*sin(psi).^5.*sin(theta).^3.*E(3,4,1)+(-30).*cos(phi).^4.*cos(psi).^3.*sin(phi).*sin(psi).^2.*sin(theta).^4.*E(3,4,1)+20.*cos(phi).^4.*cos(psi).*sin(phi).*sin(psi).^4.*sin(theta).^4.*E(3,4,1)+10.*cos(phi).^5.*cos(psi).^2.*sin(psi).^3.*sin(theta).^5.*E(3,4,1)+10.*cos(phi).^2.*cos(theta).^2.*sin(phi).^3.*sin(psi).^3.*E(4,1,3)+30.*cos(phi).^3.*cos(psi).*cos(theta).^2.*sin(phi).^2.*sin(psi).^2.*sin(theta).*E(4,1,3)+30.*cos(phi).^4.*cos(psi).^2.*cos(theta).^2.*sin(phi).*sin(psi).*sin(theta).^2.*E(4,1,3)+10.*cos(phi).^5.*cos(psi).^3.*cos(theta).^2.*sin(theta).^3.*E(4,1,3)+(-20).*cos(phi).*cos(psi).*cos(theta).*sin(phi).^4.*sin(psi).^3.*E(4,2,2)+(-60).*cos(phi).^3.*cos(psi).^3.*cos(theta).*sin(phi).^2.*sin(psi).*sin(theta).^2.*E(4,2,2)+60.*cos(phi).^3.*cos(psi).*cos(theta).*sin(phi).^2.*sin(psi).^3.*sin(theta).^2.*E(4,2,2)+(-20).*cos(phi).^4.*cos(psi).^4.*cos(theta).*sin(phi).*sin(theta).^3.*E(4,2,2)+60.*cos(phi).^4.*cos(psi).^2.*cos(theta).*sin(phi).*sin(psi).^2.*sin(theta).^3.*E(4,2,2)+20.*cos(phi).^5.*cos(psi).^3.*cos(theta).*sin(psi).*sin(theta).^4.*E(4,2,2)+(-30).*cos(phi).^2.*cos(psi).^2.*sin(phi).^3.*sin(psi).^2.*sin(2.*theta).*E(4,2,2)+10.*cos(phi).^2.*sin(phi).^3.*sin(psi).^4.*sin(2.*theta).*E(4,2,2)+10.*cos(psi).^2.*sin(phi).^5.*sin(psi).^3.*E(4,3,1)+30.*cos(phi).*cos(psi).^3.*sin(phi).^4.*sin(psi).^2.*sin(theta).*E(4,3,1)+(-20).*cos(phi).*cos(psi).*sin(phi).^4.*sin(psi).^4.*sin(theta).*E(4,3,1)+30.*cos(phi).^2.*cos(psi).^4.*sin(phi).^3.*sin(psi).*sin(theta).^2.*E(4,3,1)+(-60).*cos(phi).^2.*cos(psi).^2.*sin(phi).^3.*sin(psi).^3.*sin(theta).^2.*E(4,3,1)+10.*cos(phi).^2.*sin(phi).^3.*sin(psi).^5.*sin(theta).^2.*E(4,3,1)+10.*cos(phi).^3.*cos(psi).^5.*sin(phi).^2.*sin(theta).^3.*E(4,3,1)+(-60).*cos(phi).^3.*cos(psi).^3.*sin(phi).^2.*sin(psi).^2.*sin(theta).^3.*E(4,3,1)+30.*cos(phi).^3.*cos(psi).*sin(phi).^2.*sin(psi).^4.*sin(theta).^3.*E(4,3,1)+(-20).*cos(phi).^4.*cos(psi).^4.*sin(phi).*sin(psi).*sin(theta).^4.*E(4,3,1)+30.*cos(phi).^4.*cos(psi).^2.*sin(phi).*sin(psi).^3.*sin(theta).^4.*E(4,3,1)+10.*cos(phi).^5.*cos(psi).^3.*sin(psi).^2.*sin(theta).^5.*E(4,3,1)+5.*cos(phi).*cos(theta).*sin(phi).^4.*sin(psi).^4.*E(5,1,2)+30.*cos(phi).^3.*cos(psi).^2.*cos(theta).*sin(phi).^2.*sin(psi).^2.*sin(theta).^2.*E(5,1,2)+20.*cos(phi).^4.*cos(psi).^3.*cos(theta).*sin(phi).*sin(psi).*sin(theta).^3.*E(5,1,2)+5.*cos(phi).^5.*cos(psi).^4.*cos(theta).*sin(theta).^4.*E(5,1,2)+10.*cos(phi).^2.*cos(psi).*sin(phi).^3.*sin(psi).^3.*sin(2.*theta).*E(5,1,2)+(-5).*cos(psi).*sin(phi).^5.*sin(psi).^4.*E(5,2,1)+(-20).*cos(phi).*cos(psi).^2.*sin(phi).^4.*sin(psi).^3.*sin(theta).*E(5,2,1)+5.*cos(phi).*sin(phi).^4.*sin(psi).^5.*sin(theta).*E(5,2,1)+(-30).*cos(phi).^2.*cos(psi).^3.*sin(phi).^3.*sin(psi).^2.*sin(theta).^2.*E(5,2,1)+20.*cos(phi).^2.*cos(psi).*sin(phi).^3.*sin(psi).^4.*sin(theta).^2.*E(5,2,1)+(-20).*cos(phi).^3.*cos(psi).^4.*sin(phi).^2.*sin(psi).*sin(theta).^3.*E(5,2,1)+30.*cos(phi).^3.*cos(psi).^2.*sin(phi).^2.*sin(psi).^3.*sin(theta).^3.*E(5,2,1)+(-5).*cos(phi).^4.*cos(psi).^5.*sin(phi).*sin(theta).^4.*E(5,2,1)+20.*cos(phi).^4.*cos(psi).^3.*sin(phi).*sin(psi).^2.*sin(theta).^4.*E(5,2,1)+5.*cos(phi).^5.*cos(psi).^4.*sin(psi).*sin(theta).^5.*E(5,2,1)+sin(phi).^5.*sin(psi).^5.*E(6,1,1)+5.*cos(phi).*cos(psi).*sin(phi).^4.*sin(psi).^4.*sin(theta).*E(6,1,1)+10.*cos(phi).^2.*cos(psi).^2.*sin(phi).^3.*sin(psi).^3.*sin(theta).^2.*E(6,1,1)+10.*cos(phi).^3.*cos(psi).^3.*sin(phi).^2.*sin(psi).^2.*sin(theta).^3.*E(6,1,1)+5.*cos(phi).^4.*cos(psi).^4.*sin(phi).*sin(psi).*sin(theta).^4.*E(6,1,1)+cos(phi).^5.*cos(psi).^5.*sin(theta).^5.*E(6,1,1);
end