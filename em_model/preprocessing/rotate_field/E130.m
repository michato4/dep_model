function [res] = E130(E,phi,theta,psi)
 res=cos(theta).*sin(phi).*(cos(theta).*sin(phi).*(cos(theta).*sin(phi).*((-1).*sin(theta).*E(1,1,5)+cos(theta).*(sin(psi).*E(1,2,4)+cos(psi).*E(2,1,4)))+(cos(phi).*cos(psi)+sin(phi).*sin(psi).*sin(theta)).*((-1).*sin(theta).*E(1,2,4)+cos(theta).*(sin(psi).*E(1,3,3)+cos(psi).*E(2,2,3)))+((-1).*cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin(theta)).*((-1).*sin(theta).*E(2,1,4)+cos(theta).*(sin(psi).*E(2,2,3)+cos(psi).*E(3,1,3))))+(cos(phi).*cos(psi)+sin(phi).*sin(psi).*sin(theta)).*(cos(theta).*sin(phi).*((-1).*sin(theta).*E(1,2,4)+cos(theta).*(sin(psi).*E(1,3,3)+cos(psi).*E(2,2,3)))+(cos(phi).*cos(psi)+sin(phi).*sin(psi).*sin(theta)).*((-1).*sin(theta).*E(1,3,3)+cos(theta).*(sin(psi).*E(1,4,2)+cos(psi).*E(2,3,2)))+((-1).*cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin(theta)).*((-1).*sin(theta).*E(2,2,3)+cos(theta).*(sin(psi).*E(2,3,2)+cos(psi).*E(3,2,2))))+((-1).*cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin(theta)).*(cos(theta).*sin(phi).*((-1).*sin(theta).*E(2,1,4)+cos(theta).*(sin(psi).*E(2,2,3)+cos(psi).*E(3,1,3)))+(cos(phi).*cos(psi)+sin(phi).*sin(psi).*sin(theta)).*((-1).*sin(theta).*E(2,2,3)+cos(theta).*(sin(psi).*E(2,3,2)+cos(psi).*E(3,2,2)))+((-1).*cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin(theta)).*((-1).*sin(theta).*E(3,1,3)+cos(theta).*(sin(psi).*E(3,2,2)+cos(psi).*E(4,1,2)))))+(cos(phi).*cos(psi)+sin(phi).*sin(psi).*sin(theta)).*(cos(theta).*sin(phi).*(cos(theta).*sin(phi).*((-1).*sin(theta).*E(1,2,4)+cos(theta).*(sin(psi).*E(1,3,3)+cos(psi).*E(2,2,3)))+(cos(phi).*cos(psi)+sin(phi).*sin(psi).*sin(theta)).*((-1).*sin(theta).*E(1,3,3)+cos(theta).*(sin(psi).*E(1,4,2)+cos(psi).*E(2,3,2)))+((-1).*cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin(theta)).*((-1).*sin(theta).*E(2,2,3)+cos(theta).*(sin(psi).*E(2,3,2)+cos(psi).*E(3,2,2))))+(cos(phi).*cos(psi)+sin(phi).*sin(psi).*sin(theta)).*(cos(theta).*sin(phi).*((-1).*sin(theta).*E(1,3,3)+cos(theta).*(sin(psi).*E(1,4,2)+cos(psi).*E(2,3,2)))+(cos(phi).*cos(psi)+sin(phi).*sin(psi).*sin(theta)).*((-1).*sin(theta).*E(1,4,2)+cos(theta).*(sin(psi).*E(1,5,1)+cos(psi).*E(2,4,1)))+((-1).*cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin(theta)).*((-1).*sin(theta).*E(2,3,2)+cos(theta).*(sin(psi).*E(2,4,1)+cos(psi).*E(3,3,1))))+((-1).*cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin(theta)).*(cos(theta).*sin(phi).*((-1).*sin(theta).*E(2,2,3)+cos(theta).*(sin(psi).*E(2,3,2)+cos(psi).*E(3,2,2)))+(cos(phi).*cos(psi)+sin(phi).*sin(psi).*sin(theta)).*((-1).*sin(theta).*E(2,3,2)+cos(theta).*(sin(psi).*E(2,4,1)+cos(psi).*E(3,3,1)))+((-1).*cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin(theta)).*((-1).*sin(theta).*E(3,2,2)+cos(theta).*(sin(psi).*E(3,3,1)+cos(psi).*E(4,2,1)))))+((-1).*cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin(theta)).*(cos(theta).*sin(phi).*(cos(theta).*sin(phi).*((-1).*sin(theta).*E(2,1,4)+cos(theta).*(sin(psi).*E(2,2,3)+cos(psi).*E(3,1,3)))+(cos(phi).*cos(psi)+sin(phi).*sin(psi).*sin(theta)).*((-1).*sin(theta).*E(2,2,3)+cos(theta).*(sin(psi).*E(2,3,2)+cos(psi).*E(3,2,2)))+((-1).*cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin(theta)).*((-1).*sin(theta).*E(3,1,3)+cos(theta).*(sin(psi).*E(3,2,2)+cos(psi).*E(4,1,2))))+(cos(phi).*cos(psi)+sin(phi).*sin(psi).*sin(theta)).*(cos(theta).*sin(phi).*((-1).*sin(theta).*E(2,2,3)+cos(theta).*(sin(psi).*E(2,3,2)+cos(psi).*E(3,2,2)))+(cos(phi).*cos(psi)+sin(phi).*sin(psi).*sin(theta)).*((-1).*sin(theta).*E(2,3,2)+cos(theta).*(sin(psi).*E(2,4,1)+cos(psi).*E(3,3,1)))+((-1).*cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin(theta)).*((-1).*sin(theta).*E(3,2,2)+cos(theta).*(sin(psi).*E(3,3,1)+cos(psi).*E(4,2,1))))+((-1).*cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin(theta)).*(cos(theta).*sin(phi).*((-1).*sin(theta).*E(3,1,3)+cos(theta).*(sin(psi).*E(3,2,2)+cos(psi).*E(4,1,2)))+(cos(phi).*cos(psi)+sin(phi).*sin(psi).*sin(theta)).*((-1).*sin(theta).*E(3,2,2)+cos(theta).*(sin(psi).*E(3,3,1)+cos(psi).*E(4,2,1)))+((-1).*cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin(theta)).*((-1).*sin(theta).*E(4,1,2)+cos(theta).*(sin(psi).*E(4,2,1)+cos(psi).*E(5,1,1)))));
end