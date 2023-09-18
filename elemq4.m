% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: elemq4
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
function [EE,NDOUBGS]= elemq4(COOR,ELAS,POISN,t,PLA,WGT,NGAUSS,IPLSTR,IB,IDQ4)
NDOUBGS= NGAUSS^2;
EE= zeros(12,12);
C= zeros(3,3);
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

for I= 1:NGAUSS
    for J= 1:NGAUSS
      RI= PLA(I);
      SI= PLA(J);
      [B,VJACOB]= bmatq4(IB,COOR,IDQ4,RI,SI);
      BC= zeros(8,3);
      for K= 1:8
        for L= 1:3
          for M= 1:3
            BC(K,L)=BC(K,L)+B(M,K)*C(M,L);
          end
        end
      end
      for K= 1:8
        for L= 1:8
          for M= 1:3
            EE(K,L)=EE(K,L)+t*WGT(I)*WGT(J)*VJACOB*BC(K,M)*B(M,L);
          end
        end
      end
    end
end

end
