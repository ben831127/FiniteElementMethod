% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: stiffbc
% -----------------------------------------------------------------------------------------
% Input variables
%   NBC             : number of beam-column elements including truss elements
%   COOR            : nodal coordinates
%   VECTY           : direction of the local y-axis of beam-cloumn elements
%   IDBC            : identification matrix for beam-column elements
%   PROP            : material properties
%   SECT            : section properties
%   NADD            : address counter of each row of [K](Ji)
%   LMBC            : element location matrices of beam-column elements including truss
%   SS              : the one-dimensional stiffness
%
% Output variables
%   SS              : the one-dimensional stiffness
% -----------------------------------------------------------------------------------------
function [SS]= stiffbc(NBC,COOR,VECTY,IDBC,PROP,SECT,NADD,LMBC,SS)
for IB= 1:NBC
    [RL,ROT]= rotation(IB,COOR,VECTY,IDBC);
    I= IDBC(3,IB);
    J= IDBC(4,IB);
    ELAS= PROP(1,I);
    SHEAR= ELAS/(2*(1+PROP(2,I)));
    AREA= SECT(1,J);
    RIY= SECT(2,J);
    RIZ= SECT(3,J);
    CK= SECT(4,J);
    [EE]= elembc(ELAS,SHEAR,AREA,RIZ,RIY,CK,RL);
    [RSR]= trans(ROT,EE);
    ND= 12;
    [SS]= assem(ND,NADD,IB,RSR,LMBC,SS);
end
end
