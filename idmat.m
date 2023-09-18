% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: idmat
% -----------------------------------------------------------------------------------------
% Input variables
%   NNOD            : number of nodes
%   NFIX            : boundary conditions for nodes
% Output variables
%   NEQ             : number of equations
%   IDND            : DOFs numbering
% -----------------------------------------------------------------------------------------
function[NEQ,IDND]= idmat(NNOD,NFIX)

IDND= zeros(6,NNOD);
NEQ= 0;
for m= 1:NNOD
    for n= 1:6
      if NFIX(n,m)== 0
        NEQ=NEQ+1;
        IDND(n,m)=NEQ;
      else
        IDND(n,m)=0;
      end
    end
end
end
