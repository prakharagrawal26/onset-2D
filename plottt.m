%Date:25-Dec-2023
%Contour plots for y-z discretization 
%Magnetic field in x, varying in z + gravity and rotation along z, boussinesq run
%Using direct equations without curl and double curl

var=3; %%%%variable to be plotted in sequence ux,uy,uz,bx,by,bz,T from 1 to 7

critical_location=1; %Row number of critical eigen value from eig_val
%y-z contour plot
%eigvec=eig_vec(:,critical_location);

%eigvec=eig_vec_critical;
Uz=eigvec((var-1)*(ny+1)*(nz+1)+1:var*(ny+1)*(nz+1));
a=max(real(Uz));
b=min(real(Uz));
c=max(a,b);
Uz1=reshape(Uz,nz+1,ny+1);
[~,h]=contourf(Y,Z,real(Uz1),30,'LineStyle', 'none');
clim([-c,c]);
Lvls = h.LevelList;                                 % <— Add This Line
h.LevelList = Lvls(Lvls ~= 0);
colormap(bluewhitered);
%Set title
title('{$U_z$}','Interpreter','latex','FontSize',25)
%Set size and location of white rectangle within figure window
set(gca,'position',[.15,0.2,0.75,.65]) %1. higher value,figure moves right %2. higher value, figure moves upward %3. width of figure %4. height of figure

%X-label, X-limit, X-ticks
xl=xlabel('$y$', 'FontName', 'Times New Roman','FontSize',25,'Color','k', 'Interpreter', 'LaTeX','Rotation',0);  %'HorizontalAlignment','center','VerticalAlignment','top'
%set(xl, 'Units', 'Normalized', 'Position', [-0.22, 0.5, 0]);    %Need to modify position as required, else use above comment for alignment

%Y-label, Y-limit, Y-ticks
yl=ylabel('$z$', 'FontName', 'Times New Roman','FontSize',25,'Color','k', 'Interpreter', 'LaTeX','Rotation',0);  %'HorizontalAlignment','right','VerticalAlignment','middle'
%set(yl, 'Units', 'Normalized', 'Position', [-0.22, 0.5, 0]);    

% %x-z contour plot
% figure
% y_location=40;
% eigvec=eig_vec(:,3);
% Uz=eigvec((var-1)*(ny+1)*(nz+1)+1:var*(ny+1)*(nz+1));
% x=linspace(0,0.2,200);
% [X1,Z1]=meshgrid(x,zz);
% Uz1=reshape(Uz,nz+1,ny+1);
% Uz=Uz1(:,y_location);
% 
% for ttt=1:length(x)  %This loop copies Uz values in x direction as x variation will only come from exp(i.kx.x)
% Uz2(:,ttt)=Uz;
% end
% 
% f=exp(1j*kx*X1);
% Uz3=real((Uz2).*(f));
% [c,h]=contourf(X1,Z1,real(Uz3),20,'LineStyle', 'none');
% Lvls = h.LevelList;                                 % <— Add This Line
% h.LevelList = Lvls(Lvls ~= 0);
% colormap(bluewhitered);
% %Set title
% title('{$U_z$}','Interpreter','latex','FontSize',25)
% %Set size and location of white rectangle within figure window
% set(gca,'position',[.15,0.2,0.75,.65]) %1. higher value,figure moves right %2. higher value, figure moves upward %3. width of figure %4. height of figure
% 
% %X-label, X-limit, X-ticks
% xl=xlabel('$x$', 'FontName', 'Times New Roman','FontSize',25,'Color','k', 'Interpreter', 'LaTeX','Rotation',0);  %'HorizontalAlignment','center','VerticalAlignment','top'
% %set(xl, 'Units', 'Normalized', 'Position', [-0.22, 0.5, 0]);    %Need to modify position as required, else use above comment for alignment
% 
% %Y-label, Y-limit, Y-ticks
% yl=ylabel('$z$', 'FontName', 'Times New Roman','FontSize',25,'Color','k', 'Interpreter', 'LaTeX','Rotation',0);  %'HorizontalAlignment','right','VerticalAlignment','middle'
% %set(yl, 'Units', 'Normalized', 'Position', [-0.22, 0.5, 0]);    