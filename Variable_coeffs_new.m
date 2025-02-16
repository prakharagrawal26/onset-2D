function [c1,c2,c3,c4,c5,c6]=Variable_coeffs_new(delta,theta,m,Y,Z)

Bz=Z.*exp(-Z.^2/delta^2);
Bz=Bz./max(max(Bz));

Bzp=exp(-Z.^2/delta^2).*(1-2.*Z.^2/delta^2);
Bzp=Bzp./max(max(Bz));

rho0=(1+theta.*Y).^m;

c1=(1+theta.*Y).^-1;
c2=(1+theta.*Y).^-2;
c3=Bz./rho0;
c4=Bzp./rho0;
c5=Bz;
c6=Bzp;

c1=sparse(diag(c1(:)));
c2=sparse(diag(c2(:)));
c3=sparse(diag(c3(:)));
c4=sparse(diag(c4(:)));
c5=sparse(diag(c5(:)));
c6=sparse(diag(c6(:)));

