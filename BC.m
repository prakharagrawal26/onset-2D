function [BCUx,BCUy,BCUz,BCBx,BCBy,BCBz,BCS]=BC(ny,nz,Dy,Dz,BCzvel,BCzmag,BCyvel,BCymag,Iy2,Iz)

% Z-boundary
if BCzvel==2 && BCzmag==1   %%Stress-free and insulating with isentropy
    tt=[2,2,1,1,1,2,1];
elseif BCzvel==1 && BCzmag==1   %%No-slip and insulating with isentropy
    tt=[1,1,1,1,1,2,1];
end
BCzUx=call_BC(nz,tt(1),Dz);
BCzUy=call_BC(nz,tt(2),Dz);
BCzUz=call_BC(nz,tt(3),Dz);
BCzBx=call_BC(nz,tt(4),Dz);
BCzBy=call_BC(nz,tt(5),Dz);
BCzBz=call_BC(nz,tt(6),Dz);
BCzS=call_BC(nz,tt(7),Dz);

% Y-boundary
if BCyvel==2 && BCymag==1   %%Stress-free and insulating with isentropy
    tt=[2,1,2,1,2,1,1];
elseif BCyvel==1 && BCymag==1   %%No-slip and insulating with isentropy
    tt=[1,1,1,1,2,1,1];
end
BCyUx=call_BC(ny,tt(1),Dy);
BCyUy=call_BC(ny,tt(2),Dy);
BCyUz=call_BC(ny,tt(3),Dy);
BCyBx=call_BC(ny,tt(4),Dy);
BCyBy=call_BC(ny,tt(5),Dy);
BCyBz=call_BC(ny,tt(6),Dy);
BCyS=call_BC(ny,tt(7),Dy);

%Combine B.C.
BCzUx=kron(Iy2,sparse(BCzUx));
BCzUy=kron(Iy2,sparse(BCzUy));
BCzUz=kron(Iy2,sparse(BCzUz));
BCzBx=kron(Iy2,sparse(BCzBx));
BCzBy=kron(Iy2,sparse(BCzBy));
BCzBz=kron(Iy2,sparse(BCzBz));
BCzS=kron(Iy2,sparse(BCzS));

BCyUx=kron(sparse(BCyUx),Iz);
BCyUy=kron(sparse(BCyUy),Iz);
BCyUz=kron(sparse(BCyUz),Iz);
BCyBx=kron(sparse(BCyBx),Iz);
BCyBy=kron(sparse(BCyBy),Iz);
BCyBz=kron(sparse(BCyBz),Iz);
BCyS=kron(sparse(BCyS),Iz);

BCUx=BCzUx+BCyUx;
BCUy=BCzUy+BCyUy;
BCUz=BCzUz+BCyUz;
BCBx=BCzBx+BCyBx;
BCBy=BCzBy+BCyBy;
BCBz=BCzBz+BCyBz;
BCS=BCzS+BCyS;
end