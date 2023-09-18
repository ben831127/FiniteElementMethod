% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: load
% -----------------------------------------------------------------------------------------
% Input variables
%   NEQ             : number of equations
%   NNOD            : number of nodes
%   IDND            : DOFs numbering
%   EXLD            : external load
% Output variables
%   GLOAD           : global load
% -----------------------------------------------------------------------------------------
function[GLOAD]= load(NEQ,NNOD,IDND,EXLD)

GLOAD= zeros(NEQ,1);

for i= 1:NNOD
    for j= 1:6
      ID= IDND(j,i);
      if(ID>0)
        GLOAD(ID)= GLOAD(ID)+EXLD(j,i);
      end
    end
end

end
