% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: stifft3
% -----------------------------------------------------------------------------------------
% Input variables
%   NT3          : number of constant strain triangular elements (T3)
%   COOR         : nodal coordinates
%   IDT3         : identification matrix for T3 elements
%   PROP         : material properties
%   SECT         : section properties
%   LMT3         : element location matrices of constant strain triangular elements
%   NADD         : address counter of each row of [K](Ji)
%   IPLSTR       : indicator for plane stress ( = 1) or plane strain ( = 2)
%   SS           : the one-dimensional stiffness
%
% Output variables
%   SS           : the one-dimensional stiffness
% -----------------------------------------------------------------------------------------
function [SS]= stifft3(NT3,COOR,IDT3,PROP,SECT,LMT3,NADD,IPLSTR,SS)
for IB= 1:NT3
    I= IDT3(4,IB);
    J= IDT3(5,IB);
    ELAS= PROP(1,I);
    POISN= PROP(2,I);
    t= SECT(1,J);
    [EE]= elemt3(IB,ELAS,IDT3,COOR,POISN,t,IPLSTR);
    ND= 6;
    [SS]= assem(ND,NADD,IB,EE,LMT3,SS);
end
end