%Date-26-Dec-2023
%This code uses equations directly without taking curl and double curl
%chebyshev collocation in both y and z
%Make changes to parameters function for different parameters and B.C., no change reqd here

%Pressure has been eliminated using y-momentum equation, hence 7
%eqn,7variable, doubtful if no. of B.C. are correct, but giving correct
%results with this.

clear;
clc;
tic

[Ek,Pr,Pm,q,ny,nz,els,k1,p,Bfield,BCzmag,BCzvel,BCymag,BCyvel,Asp,delta,sigma1,Ra,m,chi]=parameters_singleRa;
theta=-chi.^(-1/m).*(-1 + chi.^(1./m));
%Identity matrices
Iy=sparse(eye(ny+1,ny+1));
Iz=sparse(eye(nz+1,nz+1));
I=kron(Iy,Iz);
Iy2=Iy; Iy2(1,:)=0; Iy2(end,:)=0;
Iz2=Iz; Iz2(1,:)=0; Iz2(end,:)=0;
I2=kron(Iy2,Iz2);

% Y-differentiation matrices
[Dy1d,yy] = cheb(ny);%multiplied with 2 to accomodate change of variable in z, compulsory else incorrect eigenfunctions
Dy1d=sparse(2*Dy1d/Asp);
Dy=Dy1d;
D2y=Dy*Dy;
yy=sparse(Asp*(yy+1)/2);
Dy=kron(Dy,Iz);
D2y=kron(D2y,Iz);

% Z-differentiation matrices
[Dz1d,zz]=cheb(nz);%multiplied with 2 to accomodate change of variable in z, compulsory else incorrect eigenfunctions
Dz1d=sparse(2*Dz1d);
Dz=Dz1d;
D2z=Dz*Dz;
Dz=kron(Iy,Dz);
D2z=kron(Iy,D2z);
zz=sparse((zz+1)/2);

[Y,Z]=meshgrid(yy,zz);

%Call B.C.
[BCUx,BCUy,BCUz,BCBx,BCBy,BCBz,BCS]=BC(ny,nz,Dy1d,Dz1d,BCzvel,BCzmag,BCyvel,BCymag,Iy2,Iz);
%Call Variable coefficients
[c1,c2,c3,c4,c5,c6]=Variable_coeffs_new(delta,theta,m,Y,Z);
eig_val=zeros(p,length(k1));
eig_vec=zeros(7*(ny+1)*(nz+1),p,length(k1));
for i=1:length(k1)
    kx=k1(i);
    disp(kx)

    %Call matrix building function
    [A,B]=build_matrix(kx,Ek,Pr,Pm,els,Ra,m,theta,BCUx,BCUy,BCUz,BCBx,BCBy,BCBz,BCS,c1,c2,c3,c4,c5,c6,I,I2,Dy,D2y,Dz,D2z);
    opts = struct('disp', 0, 'maxit', 500, 'UseParallel', true);
    [vec, D] = eigs(A,B,p,sigma1,opts);
    eig_val(:,i)=diag(D);
    eig_vec(:,:,i)=vec;
    %eig_val_sorted=sort(eig_val,'ComparisonMethod','real');
end
toc