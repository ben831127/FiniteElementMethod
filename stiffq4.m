% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: stiffq4
% -----------------------------------------------------------------------------------------
% Input variables
%   ?            : ?
%   ?            : ?
%   ?            : ?
%   ¡K
% Output variables
%   ?            : ?
% -----------------------------------------------------------------------------------------
function [SS,NDOUBGS]= stiffq4(NQ4,COOR,SS,IDQ4,PROP,SECT,PLA,WGT,NGAUSS,IPLSTR,NADD,LMQ4)
for IB= 1:NQ4
    ELAS= PROP(1,IDQ4(5,IB));
    POISN= PROP(2,IDQ4(5,IB));
    t= SECT(1,IDQ4(6,IB)); 
    [EE,NDOUBGS]= elemq4(COOR,ELAS,POISN,t,PLA,WGT,NGAUSS,IPLSTR,IB,IDQ4);
    ND= 8;
    [SS]= assem(ND,NADD,IB,EE,LMQ4,SS);
end
end
