% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: assem
% -----------------------------------------------------------------------------------------
% Input variables
%   ND              : number of displacement (12)
%   NADD            : address counter of each row of [K](Ji)
%   IB              : The ordinal number of the IB-th element              
%   RSR             : transformed stiffness matrix              
%   LMBC            : element location matrices of beam-column elements including truss             
%   SS              : the one-dimensional stiffness
% Output variables
%   SS              : the one-dimensional stiffness
% -----------------------------------------------------------------------------------------
function [SS]= assem(ND,NADD,IB,RSR,LMBC,SS)
%
% Note that we have to assemble the transformed stiffness matrix (RSR) of the IB-th 3D BC
% element into the stiffness matrix (SS, size: MSS by 1) of the whole structure using the
% location matrix of all the 3D BC elements (LMBC, size: 12 by NBC) and the address counter
% vector of each row of SS (NADD, size: NEQ by 1).
%
for M= 1:ND
    for N= 1:M
      I= LMBC(M,IB);
      J= LMBC(N,IB);
      if(I> 0 && J> 0)
        if(I>= J)
          SS(NADD(I)+J)=SS(NADD(I)+J)+RSR(M,N);
        else
          SS(NADD(J)+I)=SS(NADD(J)+I)+RSR(M,N);
        end
      end
    end
end
end
