function [Ek,Pr,Pm,q,ny,nz,els,kx,p,Bfield,BCzmag,BCzvel,BCymag,BCyvel,Asp,delta,sigma1,Racmin,m,chi]=parameters_singleRa
%Parameters
Ek=5e-5;
Pr=1;
Pm=1;
q=Pm/Pr;
ny =40;
nz =45;
%elsm=[0,0.001,0.005,0.1,0.2,0.4,1];
els=0.001;
chi=1;
%fname = sprintf('runs/chi%d_els%d.mat',chi,els);
%load(fname,'kx','Racmin');
kx=16.7894736842105;
Racmin=5.005;
p=20;  %Number of eigen values asked
BCzmag=1; %1-insulating, 2-mixed(bottom conducting,top insulating)
BCzvel=2; %1-no slip, 2-stressfree
BCymag=1; %1-insulating, 2-mixed(bottom conducting,top insulating)
BCyvel=2; %1-no slip, 2-stressfree
Bfield=2; %1-homogeneous(constant),2-inhomogeneous
Asp=1;
delta=0.14;
sigma1=1e-6;
m=1.495;
end