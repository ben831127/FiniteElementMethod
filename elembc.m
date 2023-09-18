% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: elembc
% -----------------------------------------------------------------------------------------
% Input variables
%   ELAS           : Young's modulus
%   SHEAR          : SHEAR modulus
%   AREA           : area of the the element section
%   RIZ            : area moment of inertia of Z-axia
%   RIY            : area moment of inertia of Y-axia
%   CK             : second axial moment of area
%   RL             : Length of the IB-th 3D BC element
% Output variables
%   EE             : local stiffness matrix
% -----------------------------------------------------------------------------------------
function [EE]= elembc(ELAS,SHEAR,AREA,RIZ,RIY,CK,RL)
EE= zeros(12,12);
  RLSQ= RL*RL;
  RLCU= RLSQ*RL;
  GKL= SHEAR*CK/RL;
  
  EE(1,1)   =ELAS*AREA/RL;
  EE(2,2)   =12*ELAS*RIZ/RLCU;
  EE(3,3)   =12*ELAS*RIY/RLCU;
  EE(4,4)   =GKL;
  EE(3,5)   =-6*ELAS*RIY/RLSQ;
  EE(5,5)   =4*ELAS*RIY/RL;
  EE(2,6)   =6*ELAS*RIZ/RLSQ;
  EE(6,6)   =4*ELAS*RIZ/RL;
  EE(1,7)   =-EE(1,1);
  EE(7,7)   =EE(1,1);
  EE(2,8)   =-EE(2,2);
  EE(6,8)   =-EE(2,6);
  EE(8,8)   =EE(2,2);
  EE(3,9)   =-EE(3,3);
  EE(5,9)   =-EE(3,5);
  EE(9,9)   =EE(3,3);
  EE(4,10)  =-EE(4,4);
  EE(10,10) =EE(4,4);
  EE(3,11)  =EE(3,5);
  EE(5,11)  =EE(5,5)/2;
  EE(9,11)  =EE(5,9);
  EE(11,11) =EE(5,5);
  EE(2,12)  =EE(2,6);
  EE(6,12)  =EE(6,6)/2;
  EE(8,12)  =EE(6,8);
  EE(12,12) =EE(6,6);
  
  for I=1:11
    for J=I+1:12
	  EE(J,I)=EE(I,J);
    end
  end
end
