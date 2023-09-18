% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: rotation
% -----------------------------------------------------------------------------------------
% Input variables
%   IB           : The ordinal number of the IB-th element
%   COOR         : nodal coordinates
%   VECTY        : direction of the local y-axis of beam-cloumn elements
%   IDBC         : identification matrix for beam-column elements
% Output variables
%   RL           : Length of the IB-th 3D BC element
%   ROT          : Rotation matrix of the IB-th 3D BC element
% -----------------------------------------------------------------------------------------
function [RL,ROT]= rotation(IB,COOR,VECTY,IDBC)
ROT= zeros(3,3);
  for I= 1:3
    ROT(1,I)= COOR(I,IDBC(2,IB))-COOR(I,IDBC(1,IB));
    ROT(2,I)= VECTY(I,IB);
  end
  ROT(3,1)= ROT(1,2)*ROT(2,3)-ROT(1,3)*ROT(2,2);
  ROT(3,2)= ROT(1,3)*ROT(2,1)-ROT(1,1)*ROT(2,3);
  ROT(3,3)= ROT(1,1)*ROT(2,2)-ROT(1,2)*ROT(2,1);
%
% Note the case that some of row vectors in the rotation matrix may be the (0,0,0) vector,
% which length is zero, while dealing with the corresponding unit vector.
%
  L1=0;
  L2=0;
  L3=0;
  for I=1:3
    L1=ROT(1,I)^2+L1;
    L2=ROT(2,I)^2+L2;
	L3=ROT(3,I)^2+L3;
  end
      
  LX=L1^0.5;
  LY=L2^0.5;
  LZ=L3^0.5;
  RL= LX;
  for I=1:3
	ROT(1,I)=ROT(1,I)/LX;

    if (LY~=0)
	  ROT(2,I)=ROT(2,I)/LY;
    end
	if (LZ~=0)
	  ROT(3,I)=ROT(3,I)/LZ;
	end
  end
end