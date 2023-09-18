% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: forcet3
% -----------------------------------------------------------------------------------------
% Input variables
%   IPLSTR      : ?
%   ?            : ?
%   ?            : ?
%   ¡K
% Output variables
%   ?            : ?
% -----------------------------------------------------------------------------------------
function [STRAINT3,STRESST3]= forcet3(IPLSTR,NT3,COOR,IDT3,PROP,LMT3,DELTA)
  STRAINT3= zeros(3,NT3);
  STRESST3= zeros(3,NT3);
  
  for IB= 1:NT3
    DSG= zeros(6);
    for I= 1:6
      if(LMT3(I,IB) > 0)
        DSG(I)= DELTA(LMT3(I,IB));
      else 
        DSG(I)= 0;
      end 
    end

    X= zeros(3);
    Y= zeros(3);
    for I= 1:3
      X(I)= COOR(1,IDT3(I,IB));
      Y(I)= COOR(2,IDT3(I,IB));
    end

    AE= 0.5*(X(2)*Y(3)-X(3)*Y(2)+X(3)*Y(1)-X(1)*Y(3)+X(1)*Y(2)-X(2)*Y(1));
    B(1,1)= (Y(2)-Y(3))/(2*AE); 
    B(1,2)= 0;
    B(1,3)= (Y(3)-Y(1))/(2*AE) ;
    B(1,4)= 0;
    B(1,5)=(Y(1)-Y(2))/(2*AE) ;
    B(1,6)=0;
    B(2,1)=0;
    B(2,2)=(X(3)-X(2))/(2*AE);
    B(2,3)=0;
    B(2,4)=(X(1)-X(3))/(2*AE);
    B(2,5)=0;
    B(2,6)=(X(2)-X(1))/(2*AE);
    B(3,1)=B(2,2);
    B(3,2)=B(1,1);
    B(3,3)=B(2,4);
    B(3,4)=B(1,3);
    B(3,5)=B(2,6);
    B(3,6)=B(1,5);

    ELAS= PROP(1,IDT3(4,IB));
    NU= PROP(2,IDT3(4,IB));
    if(IPLSTR== 1)
      C(1,1)= ELAS/(1-NU^2);
      C(1,2)= C(1,1)*NU;
      C(1,3)= 0;
      C(2,2)= C(1,1);
      C(2,3)= 0;
      C(3,3)= C(1,1)*(1-NU)/2;
    else
      C(1,1)= ELAS*(1-NU)/(1+NU)/(1-2*NU);
      C(1,2)= C(1,1)*NU/(1-NU);
      C(1,3)= 0;
      C(2,2)= C(1,1);
      C(2,3)= 0;
      C(3,3)= ELAS/2/(1+NU);
    end
    C(2,1)= C(1,2);
    C(3,1)= C(1,3);
    C(3,2)= C(2,3);

    for I= 1:3
      STRAINT3(I,IB)= 0;
      for J= 1:6
        STRAINT3(I,IB)=STRAINT3(I,IB)+B(I,J)*DSG(J);
      end
    end

    for I= 1:3
      STRESST3(I,IB)= 0;
      for J= 1:3
        STRESST3(I,IB)=STRESST3(I,IB)+C(I,J)*STRAINT3(J,IB);
      end
    end
  end
end
