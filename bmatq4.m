% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: bmatq4
% -----------------------------------------------------------------------------------------
function[B,VJACOB]= bmatq4(IB,COOR,IDQ4,RI,SI)
  JACOB= zeros(2,2);
JACOB(1,1)= ...
0.25*((COOR(1,IDQ4(2,IB))-COOR(1,IDQ4(1,IB)))*(1-SI)+...
(COOR(1,IDQ4(3,IB))-COOR(1,IDQ4(4,IB)))*(1+SI));
  JACOB(1,2)= ...
0.25*((COOR(2,IDQ4(2,IB))-COOR(2,IDQ4(1,IB)))*(1-SI)+...
(COOR(2,IDQ4(3,IB))-COOR(2,IDQ4(4,IB)))*(1+SI));
JACOB(2,1)= ...
0.25*((COOR(1,IDQ4(3,IB))-COOR(1,IDQ4(2,IB)))*(1+RI)+...
(COOR(1,IDQ4(4,IB))-COOR(1,IDQ4(1,IB)))*(1-RI));
JACOB(2,2)= ...
0.25*((COOR(2,IDQ4(3,IB))-COOR(2,IDQ4(2,IB)))*(1+RI)+...
(COOR(2,IDQ4(4,IB))-COOR(2,IDQ4(1,IB)))*(1-RI));
  VJACOB= JACOB(1,1)*JACOB(2,2)-JACOB(1,2)*JACOB(2,1);

JACOBIV= zeros(2,2);
  JACOBIV(1,1)= JACOB(2,2)/VJACOB;
  JACOBIV(1,2)= (-1)*JACOBIV(1,2)/VJACOB;
  JACOBIV(2,1)= (-1)*JACOBIV(2,1)/VJACOB;
  JACOBIV(2,2)= JACOB(1,1)/VJACOB;

DRSHAPE= zeros(2,4);
  DRSHAPE(1,1)= (-1)*0.25*(1-SI);
  DRSHAPE(1,2)= (-1)*DRSHAPE(1,1);
  DRSHAPE(1,3)= 0.25*(1+SI);
  DRSHAPE(1,4)= (-1)*DRSHAPE(1,3);
  DRSHAPE(2,1)= (-1)*0.25*(1-RI);
  DRSHAPE(2,2)= (-1)*0.25*(1+RI);
  DRSHAPE(2,3)= (-1)*DRSHAPE(2,2);
  DRSHAPE(2,4)= (-1)*DRSHAPE(2,1);

  DXSHAPE= zeros(2,4);
  for I=1:2
    for J=1:4
      for K=1:2
        DXSHAPE(I,J)=DXSHAPE(I,J)+JACOBIV(I,K)*DRSHAPE(K,J);
      end
    end
  end

  B= zeros(3,8);
  B(1,1)= DXSHAPE(1,1);
  B(1,2)= 0;
  B(1,3)= DXSHAPE(1,2);
  B(1,4)= 0;
  B(1,5)= DXSHAPE(1,3);
  B(1,6)= 0;
  B(1,7)= DXSHAPE(1,4);
  B(1,8)= 0;
  B(2,1)= 0;
  B(2,2)= DXSHAPE(2,1);
  B(2,3)= 0;
  B(2,4)= DXSHAPE(2,2);
  B(2,5)= 0;
  B(2,6)= DXSHAPE(2,3);
  B(2,7)= 0;
  B(2,8)= DXSHAPE(2,4);
  B(3,1)=DXSHAPE(2,1);
  B(3,2)=DXSHAPE(1,1);
  B(3,3)=DXSHAPE(2,2);
  B(3,4)=DXSHAPE(1,2);
  B(3,5)=DXSHAPE(2,3);
  B(3,6)=DXSHAPE(1,3);
  B(3,7)=DXSHAPE(2,4);
  B(3,8)=DXSHAPE(1,4);

end
