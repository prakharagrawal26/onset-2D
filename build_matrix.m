function [A,B]=build_matrix(kx,Ek,Pr,Pm,els,Ra,m,theta,BCUx,BCUy,BCUz,BCBx,BCBy,BCBz,BCS,c1,c2,c3,c4,c5,c6,I,I2,Dy,D2y,Dz,D2z)

D2=D2y+D2z-I*kx^2;
% %x-momentum equation: substituted pressure from this equation, hence this eqn is reduced
% 
% %y-momentum equation:
% A2=[I2*(-I+Ek*1j/kx*Dy*D2+Ek*1j*m*theta/kx*Dy*c1*Dy) I2*(Ek*D2+5/3*Ek*m*theta*c1*Dy+1j/kx*Dy+Ek*(2*m+1)*m*theta^2/3*c2-2*Ek*m*theta/3*Dy*c1)+BCUy I2*(0*I) -I2*els*c3*Dy I2*els*c3*1j*kx I2*(els*1j/kx*Dy*c4) I2*Pm/Pr*Ra];
% B2=[I2*(Ek/Pm*1j/kx*Dy) I2*(Ek/Pm) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I)];
% 
% %z-momentum equation
% A3=[I2*(1j/kx*Ek*Dz*D2+1j/kx*Ek*m*theta*Dz*c1*Dy) I2*(1j/kx*Dz) I2*(Ek*D2+Ek*m*theta*c1*Dy)+BCUz I2*(-els*c3*Dz+els*c4) I2*(0*I) I2*(1j*kx*els*c3-els*1j/kx*Dz*c4) I2*(0*I)];
% B3=[I2*(Ek*1j/(kx*Pm)*Dz) I2*(0*I) I2*(Ek/Pm*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I)];
% 
% %x-induction equation
% A4=[I2*(0*I) I2*(-c5*Dy) I2*(-c5*Dz-c6) I2*(D2)+BCBx I2*(0*I) I2*(0*I) I2*(0*I)];
% B4=[I2*(0*I) I2*(0*I) I2*(0*I) I2*(I) I2*(0*I) I2*(0*I) I2*(0*I)];
% 
% %y-induction equation
% A5=[I2*(0*I) I2*(1j*kx*c5) I2*(0*I) I2*(0*I) I2*(D2)+BCBy I2*(0*I) I2*(0*I)];
% B5=[I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(I) I2*(0*I) I2*(0*I)];
% 
% %z-induction equation
% A6=[I2*(0*I) I2*(0*I) I2*(1j*kx*c5) I2*(0*I) I2*(0*I) I2*(D2)+BCBz I2*(0*I)];
% B6=[I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(I) I2*(0*I)];
% 
% %continuity equation
% A7=[I2*(1j*kx*I)+BCUx I2*(Dy+m*theta*c1) I2*(Dz) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I)];
% B7=[I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I)];
% 
% %entropy equation
% A8=[I2*(0*I) I2*(c1) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(Pm/Pr*D2+Pm/Pr*(m+1)*theta*c1*Dy)+BCS];
% B8=[I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(I)];
% 
% A=[A2;A3;A4;A5;A6;A7;A8];
% B=[B2;B3;B4;B5;B6;B7;B8];
A=[I2*(-I+Ek*1j/kx*Dy*D2+Ek*1j*m*theta/kx*Dy*c1*Dy) I2*(Ek*D2+5/3*Ek*m*theta*c1*Dy+1j/kx*Dy+Ek*(2*m+1)*m*theta^2/3*c2-2*Ek*m*theta/3*Dy*c1)+BCUy I2*(0*I) -I2*els*c3*Dy I2*els*c3*1j*kx I2*(els*1j/kx*Dy*c4) I2*Pm/Pr*Ra; 
    I2*(1j/kx*Ek*Dz*D2+1j/kx*Ek*m*theta*Dz*c1*Dy) I2*(1j/kx*Dz) I2*(Ek*D2+Ek*m*theta*c1*Dy)+BCUz I2*(-els*c3*Dz+els*c4) I2*(0*I) I2*(1j*kx*els*c3-els*1j/kx*Dz*c4) I2*(0*I);
    I2*(0*I) I2*(-c5*Dy) I2*(-c5*Dz-c6) I2*(D2)+BCBx I2*(0*I) I2*(0*I) I2*(0*I);
    I2*(0*I) I2*(1j*kx*c5) I2*(0*I) I2*(0*I) I2*(D2)+BCBy I2*(0*I) I2*(0*I);
    I2*(0*I) I2*(0*I) I2*(1j*kx*c5) I2*(0*I) I2*(0*I) I2*(D2)+BCBz I2*(0*I);
    I2*(1j*kx*I)+BCUx I2*(Dy+m*theta*c1) I2*(Dz) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I);
    I2*(0*I) I2*(c1) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(Pm/Pr*D2+Pm/Pr*(m+1)*theta*c1*Dy)+BCS];

B=[I2*(Ek/Pm*1j/kx*Dy) I2*(Ek/Pm) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I);
    I2*(Ek*1j/(kx*Pm)*Dz) I2*(0*I) I2*(Ek/Pm*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I);
    I2*(0*I) I2*(0*I) I2*(0*I) I2*(I) I2*(0*I) I2*(0*I) I2*(0*I);
    I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(I) I2*(0*I) I2*(0*I);
    I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(I) I2*(0*I);
    I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I);
    I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(0*I) I2*(I)];