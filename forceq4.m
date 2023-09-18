%-----------------------------------------------------------------------------------------
%                             SUBROUTINE: forceq4
% -----------------------------------------------------------------------------------------
% Input variables
%   IB           : The ordinal number of the IB-th element
%   ?            : ?
%   ?            : ?
%   ?            : ?
%   ¡K
% Output variables
%   ?            : ?
% -----------------------------------------------------------------------------------------
function [STRAINQ4,STRESSQ4,COORGAUS]= forceq4(NGAUSS,NDOUBGS,NQ4,IPLSTR,COOR,IDQ4,PROP,...
  LMQ4,DELTA)
  [PLA,~]= gauss(NGAUSS);
  STRAINQ4= zeros(3*NDOUBGS+3,NQ4);
  STRESSQ4= zeros(3*NDOUBGS+3,NQ4);
  COORGAUS= zeros(3,NDOUBGS*NQ4);

  for IB= 1:NQ4
    DSG= zeros(8);
    for I= 1:8
      if(LMQ4(I,IB) > 0)
        DSG(I)=DELTA(LMQ4(I,IB));
      else 
        DSG(I)=0;
      end
    end
    X= zeros(4);
    Y= zeros(4);
    Z= zeros(4);
    for I= 1:4
      X(I)= COOR(1,IDQ4(I,IB));
      Y(I)= COOR(2,IDQ4(I,IB));
      Z(I)= COOR(3,IDQ4(I,IB));
    end

  ELAS= PROP(1,IDQ4(5,IB));
  NU= PROP(2,IDQ4(5,IB));
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

I= 1;
for M= 1:NGAUSS
      for N= 1:NGAUSS
        RI= PLA(N);
        SI= PLA(M);
        [B,~]= bmatq4(IB,COOR,IDQ4,RI,SI);
        for K= 1:3
          if(K== 1)
            for J= 1:8
              STRAINQ4(I,IB)= STRAINQ4(I,IB)+B(K,J)*DSG(J);
            end
          end
          if(K== 2)
            for J= 1:8
              STRAINQ4(I+1,IB)=STRAINQ4(I+1,IB)+B(K,J)*DSG(J);
            end
          end
          if(K== 3)
            for J= 1:8
              STRAINQ4(I+2,IB)=STRAINQ4(I+2,IB)+B(K,J)*DSG(J);
            end
          end
        end
        STRAINQ4(NDOUBGS*3+1,IB)= STRAINQ4(NDOUBGS*3+1,IB)+STRAINQ4(I,IB);
	    STRAINQ4(NDOUBGS*3+2,IB)= STRAINQ4(NDOUBGS*3+2,IB)+STRAINQ4(I+1,IB);
        STRAINQ4(NDOUBGS*3+3,IB)= STRAINQ4(NDOUBGS*3+3,IB)+STRAINQ4(I+2,IB);  
        I=I+3; 
      end
end
STRAINQ4(NDOUBGS*3+1,IB)= STRAINQ4(NDOUBGS*3+1,IB)/NDOUBGS;
STRAINQ4(NDOUBGS*3+2,IB)= STRAINQ4(NDOUBGS*3+2,IB)/NDOUBGS;
STRAINQ4(NDOUBGS*3+3,IB)= STRAINQ4(NDOUBGS*3+3,IB)/NDOUBGS;

I= 1;
    for M= 1:NGAUSS
      for N= 1:NGAUSS
        for K= 1:3
          if(K== 1)
            for J= 1:3
              STRESSQ4(I,IB)= STRESSQ4(I,IB)+C(K,J)*STRAINQ4(I+(J-1),IB);
            end
          end
          if(K== 2)
            for J= 1:3
              STRESSQ4(I+1,IB)=STRESSQ4(I+1,IB)+C(K,J)*STRAINQ4(I+(J-1),IB);
            end
          end
          if(K== 3)
            for J= 1:3
              STRESSQ4(I+2,IB)=STRESSQ4(I+2,IB)+C(K,J)*STRAINQ4(I+(J-1),IB);
            end
          end
        end
        STRESSQ4(NDOUBGS*3+1,IB)= STRESSQ4(NDOUBGS*3+1,IB)+STRESSQ4(I,IB);
        STRESSQ4(NDOUBGS*3+2,IB)= STRESSQ4(NDOUBGS*3+2,IB)+STRESSQ4(I+1,IB);
        STRESSQ4(NDOUBGS*3+3,IB)= STRESSQ4(NDOUBGS*3+3,IB)+STRESSQ4(I+2,IB);
     	  I= I+3;
      end
    end
    STRESSQ4(NDOUBGS*3+1,IB)=STRESSQ4(NDOUBGS*3+1,IB)/NDOUBGS;
    STRESSQ4(NDOUBGS*3+2,IB)=STRESSQ4(NDOUBGS*3+2,IB)/NDOUBGS;
    STRESSQ4(NDOUBGS*3+3,IB)=STRESSQ4(NDOUBGS*3+3,IB)/NDOUBGS;

    I= 4*IB-3;
    for M= 1:NGAUSS
      SI= PLA(M);
      for N= 1:NGAUSS
        RI= PLA(N);
        NSHAPE(1)= 0.25*(1-RI)*(1-SI);
        NSHAPE(2)=0.25*(1+RI)*(1-SI);
        NSHAPE(3)=0.25*(1+RI)*(1+SI);
        NSHAPE(4)=0.25*(1-RI)*(1+SI);
        for J= 1:4
          COORGAUS(1,I)=COORGAUS(1,I)+X(J)*NSHAPE(J);
          COORGAUS(2,I)=COORGAUS(2,I)+Y(J)*NSHAPE(J);
          COORGAUS(3,I)=COORGAUS(3,I)+Z(J)*NSHAPE(J);
        end
        I= I+1;
      end
    end
  end
end
