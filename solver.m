% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: solver
% -----------------------------------------------------------------------------------------
function [X]= solver(SS,X,NEQ,MHB,NLC,NONZ,NADD)
for I= 2:NEQ
if(NONZ(I)< I)
    for J= NONZ(I):I-1
      for MD=1:NLC
        X(I,MD)=X(I,MD)-SS(NADD(I)+J)*X(J,MD);
      end
    end
end
end

for I= 1:NEQ
NPS=NADD(I)+I;
if(abs(SS(NPS))< 10^(-12))
    SS(NPS)= 10^(-12)*sign(SS(NPS));
end
for MD= 1:NLC
    X(I,MD)= X(I,MD)/SS(NPS);
end
end
for I= NEQ-1:-1:1
MIB= min(I+MHB-1,NEQ);
for J= I+1:MIB
    if(NONZ(J)<= I)
      for MD= 1:NLC
        X(I,MD)= X(I,MD)-SS(NADD(J)+I)*X(J,MD);
      end
    end
end
end
end
