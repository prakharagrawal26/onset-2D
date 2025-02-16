%This function is only used to build Boundary conditions. It is called only inside BC.m
function M=call_BC(n,tt,DD) %n=ny or nz, tt=number related to given B.C., DD=first derivative matrix
M=zeros(n+1);
if tt==1
    M(1,1)=1;
    M(end,end)=1;
elseif tt==2
    M(1,:)=DD(1,:);
    M(end,:)=DD(end,:);
end
end