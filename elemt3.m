% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: elemt3
% -----------------------------------------------------------------------------------------
% Input variables
%   IB           : The ordinal number of the IB-th element
%   ELAS         : Young's modulus
%   IDT3         : identification matrix for T3 elements
%   COOR         : nodal coordinates
%   POISN        : Poisson's ratio
%   t            : thickness of the element
%   IPLSTR       : indicator for plane stress ( = 1) or plane strain ( = 2)
%
% Output variables
%   EE            : local stiffness matrix
% -----------------------------------------------------------------------------------------
function [EE]= elemt3(IB,ELAS,IDT3,COOR,POISN,t,IPLSTR)
  EE= zeros(6,6);
  if(IPLSTR== 1)
    C(1,1)=ELAS/(1-POISN*POISN)*1;
    C(1,2)=ELAS/(1-POISN*POISN)*POISN;
    C(1,3)=0;
    C(2,1)=ELAS/(1-POISN*POISN)*POISN;
    C(2,2)=ELAS/(1-POISN*POISN)*1;
    C(2,3)=0;
    C(3,1)=0;
    C(3,2)=0;
    C(3,3)=ELAS/(1-POISN*POISN)*(1-POISN)/2;
  end
  if(IPLSTR== 2)
	C(1,1)=ELAS/(1-POISN*POISN)*(1-POISN);
	C(1,2)=ELAS/(1-POISN*POISN)*POISN;
	C(1,3)=0;
	C(2,1)=ELAS/(1-POISN*POISN)*POISN;
	C(2,2)=ELAS/(1-POISN*POISN)*(1-POISN);
	C(2,3)=0;
	C(3,1)=0;
	C(3,2)=0;
	C(3,3)=ELAS/(1-POISN*POISN)*(1-2*POISN)/2;
  end
  X= zeros(3,1);
  Y= zeros(3,1);
  Z= zeros(3,1);
  for I= 1:3
	X(I)=COOR(1,IDT3(I,IB));
	Y(I)=COOR(2,IDT3(I,IB));
	Z(I)=COOR(3,IDT3(I,IB));
  end
  
  Ae= 0.5*abs(X(2)*Y(3)-X(3)*Y(2)+X(3)*Y(1)-X(1)*Y(3)+X(1)*Y(2)-X(2)*Y(1));
  
  b1=Y(2)-Y(3);
  b2=Y(3)-Y(1);
  b3=Y(1)-Y(2);
  c1=X(3)-X(2);
  c2=X(1)-X(3);
  c3=X(2)-X(1);
B= zeros(3,6);
  B(1,1)= 1/(2*Ae)*b1;
  B(1,3)= 1/(2*Ae)*b2;
  B(1,5)= 1/(2*Ae)*b3;
  B(2,2)= 1/(2*Ae)*c1;
  B(2,4)= 1/(2*Ae)*c2;
  B(2,6)= 1/(2*Ae)*c3;
  B(3,1)= 1/(2*Ae)*c1;
  B(3,2)= 1/(2*Ae)*b1;
  B(3,3)= 1/(2*Ae)*c2;
  B(3,4)= 1/(2*Ae)*b2;
  B(3,5)= 1/(2*Ae)*c3;
  B(3,6)= 1/(2*Ae)*b3;

EE= t*Ae*B.'*C*B;

end
