% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: trans
% -----------------------------------------------------------------------------------------
% Input variables
%   ROT            : Rotation matrix of the IB-th 3D BC element
%   EE             : local stiffness matrix
% Output variables
%   RSR            : global stiffness matrix
% -----------------------------------------------------------------------------------------
function [RSR]= trans(ROT,EE)
%
% Note that we have to build the transformation matrix T (TRA) by using the rotation matrix
% (ROT), and then build transformed stiffness matrix k (RSR) through the formula
% k= T^T k^¡¦ T.  (Chapter 4.5; EE represents k.)
%
TRA= ROT.';
A= zeros(3,3);
RSR= zeros(12,12);
  for I=1:4
    for J=I:4
      for K=1:3
        for L=1:3
          A(K,L)=EE(K+I*3-3,L+J*3-3);
        end
      end
	  TMP=A*ROT;
	  C=TRA*TMP;
      for K=1:3
        for L=1:3
          RSR(K+I*3-3,L+J*3-3)=C(K,L);
        end
      end
    end
  end
      
  for I=1:11
    for J=I+1:12
	  RSR(J,I)=RSR(I,J);
    end
  end
end
