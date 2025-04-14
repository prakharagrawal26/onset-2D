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

[Ek,Pr,Pm,q,ny,nz,elsm,k1,p,Bfield,BCzmag,BCzvel,BCymag,BCyvel,Asp,delta,sigma1,Ra,m,chim]=parameters;

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

for i=1:length(elsm)
    for j=1:length(chim)
        els=elsm(i);
        chi=chim(j);
        theta=-chi.^(-1/m).*(-1 + chi.^(1./m));
        %Call B.C.
        [BCUx,BCUy,BCUz,BCBx,BCBy,BCBz,BCS]=BC(ny,nz,Dy1d,Dz1d,BCzvel,BCzmag,BCyvel,BCymag,Iy2,Iz);
        %Call Variable coefficients
        [c1,c2,c3,c4,c5,c6]=Variable_coeffs_new(delta,theta,m,Y,Z);
        eig_val=zeros(p,length(k1));
        eig_vec=zeros(7*(ny+1)*(nz+1),length(k1));
        parfor rr=1:length(k1)
            kx=k1(rr);
            disp(kx)

            %Call matrix building function
            [A,B]=build_matrix(kx,Ek,Pr,Pm,els,Ra,m,theta,BCUx,BCUy,BCUz,BCBx,BCBy,BCBz,BCS,c1,c2,c3,c4,c5,c6,I,I2,Dy,D2y,Dz,D2z);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Binary Search Loop
            strt=0;
            endd=2000;
            num=0;
            D2=D2y+D2z-I*kx^2;
            A2=I2*Pm/Pr*strt;
            A(1:(ny+1)*(nz+1),6*(ny+1)*(nz+1)+1:7*(ny+1)*(nz+1))=A2;
            [~,~,flag1,lm1]=solver_new(A,B,num,p,sigma1);
            A2=I2*Pm/Pr*endd;
            A(1:(ny+1)*(nz+1),6*(ny+1)*(nz+1)+1:7*(ny+1)*(nz+1))=A2;
            [~,~,flag2,lm2]=solver_new(A,B,num,p,sigma1);
            diff=endd-strt;
            fl=1;
            while diff>0.01
                if flag1==1 && flag2==1
                    fl=0;
                    break;
                elseif flag1==0 && flag2==0
                    strt=endd;
                    endd=endd+500;
                    A2=I2*Pm/Pr*endd;
                    A(1:(ny+1)*(nz+1),6*(ny+1)*(nz+1)+1:7*(ny+1)*(nz+1))=A2;
                    [~,~,flag2,lm2]=solver_new(A,B,num,p,sigma1);
                elseif flag1==0 && flag2==1
                    midd=(strt+endd)/2;
                    A2=I2*Pm/Pr*midd;
                    A(1:(ny+1)*(nz+1),6*(ny+1)*(nz+1)+1:7*(ny+1)*(nz+1))=A2;
                    [~,~,flag3,lm3]=solver_new(A,B,num,p,sigma1);
                    if flag3==0
                        strt=midd;
                        lm1=lm3;
                        flag1=flag3;
                    elseif flag3==1
                        endd=midd;
                        lm2=lm3;
                        flag2=flag3;
                    end
                    diff=endd-strt;
                else
                    fl=2;
                    break;
                end
            end
            if fl==0
                eig_val(:,rr)=1e20;
                Rac(rr)=1e20;
            elseif fl==2
                eig_val(:,rr)=-1e20;
                Rac(rr)=-1e20;
            else
               
                Rac(rr)=endd;
                A2=I2*Pm/Pr*endd;
                A(1:(ny+1)*(nz+1),6*(ny+1)*(nz+1)+1:7*(ny+1)*(nz+1))=A2;
                [~,eig_vec1,~,~,all_eig]=solver_new(A,B,1,p,sigma1);
                 eig_val(:,rr)=all_eig;
                eig_vec(:,rr)=eig_vec1(:,1);
            end
        end
        [Racmin,index]=min(Rac);
        Racmax=max(Rac);
        kx=k1(index);
        %%% Eigenvector calculation for critical value
        [A,B]=build_matrix(kx,Ek,Pr,Pm,els,Racmin,m,theta,BCUx,BCUy,BCUz,BCBx,BCBy,BCBz,BCS,c1,c2,c3,c4,c5,c6,I,I2,Dy,D2y,Dz,D2z);
        [eig_val_critical,eig_vec_critical,fl,lm]=solver_new(A,B,1,p,sigma1);
        
        clear A B
        clear BCUx BCUy BCUz BCBx BCBy BCBz BCS
        clear c1 c2 c3 c4 c5 c6

        if Racmax==1e20
            a='error found';
        end
        if Racmin==-1e20
            aa='other error found';
        end
        toc
        fname = sprintf('runs/chi%d_els%d.mat',chi,els);
        save(fname)
    end
    
end
toc
