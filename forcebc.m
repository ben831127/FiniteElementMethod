% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: forcebc
% -----------------------------------------------------------------------------------------
% Input variables
%   NBC          : ?
%   ?            : ?
%   ?            : ?
%   ¡K
% Output variables
%   ELFOR       : ?
% -----------------------------------------------------------------------------------------
function [ELFOR]= forcebc(NBC,COOR,IDBC,VECTY,PROP,SECT,LMBC,DELTA)
  ELFOR= zeros(12,NBC);
  for IB= 1:NBC
    DSG= zeros(12);
    DSL= zeros(12);
    for I= 1:12
      if(LMBC(I,IB)> 0)
        DSG(I)= DELTA(LMBC(I,IB));
      else 
        DSG(I)= 0;
      end
    end
    [RL,ROT]= ROTATION(IB,COOR,VECTY,IDBC);

    I=1; 
    for K= 1:3
      for J= 1:3
        DSL(I)= DSL(I)+ROT(K,J)*DSG(J);
        DSL(I+3)= DSL(I+3)+ROT(K,J)*DSG(J+3);
        DSL(I+6)= DSL(I+6)+ROT(K,J)*DSG(J+6);
        DSL(I+9)= DSL(I+9)+ROT(K,J)*DSG(J+9);
      end
	I=I+1;
    end

ELAS= PROP(1,IDBC(3,IB));
SHEAR= ELAS/2/(1+PROP(2,IDBC(3,IB)));
AREA= SECT(1,IDBC(4,IB));
RIY= SECT(2,IDBC(4,IB));
RIZ= SECT(3,IDBC(4,IB));
CK= SECT(4,IDBC(4,IB));
    [EE]= ELEMBC(ELAS,SHEAR,AREA,RIZ,RIY,CK,RL);
    for I= 1:12
      ELFOR(I,IB)= 0;
      for J= 1:12
        ELFOR(I,IB)= ELFOR(I,IB)+EE(I,J)*DSL(J);
      end
    end
  end
end
