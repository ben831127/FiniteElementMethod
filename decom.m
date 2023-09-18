% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: decom
% -----------------------------------------------------------------------------------------
function[S,ISIGN]= decom(NONZ,NADD,S,MHB,NEQ)
ISIGN=0;
if(S(1)== 0)
    ISIGN= 1;
    return;
end

for I= 2:MHB
    if(NONZ(I)== 1)
      I1= NADD(I)+1;
      S(I1)= S(I1)/S(1);
    end
end
for J= 2:NEQ
    JA= NADD(J)+J;
if(NONZ(J)< J)
      for L= NONZ(J):J-1
        JB= NADD(J)+L;
        S(JA)= S(JA)-S(NADD(L)+L)*S(JB)*S(JB);
      end
end
if(S(JA)== 0)
      ISIGN= 1;
      return;    
end
    if(J== NEQ) 
      return;
    end
    MIB= min(J+MHB-1,NEQ);
if(MIB>= J+1)
      for K= J+1:MIB
        if(NONZ(K)<= J)
          KA= NADD(K)+J;
          MX= max(NONZ(K),NONZ(J));
          if(MX< J)
            for L= MX:J-1
              S(KA)=S(KA)-S(NADD(L)+L)*S(NADD(J)+L)*S(NADD(K)+L);
            end
          end
          S(KA)= S(KA)/S(JA);
        end
      end
end 
end
end
