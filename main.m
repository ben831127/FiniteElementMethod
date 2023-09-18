% -----------------------------------------------------------------------------------------
%                             PROGRAM OF FINITE ELEMENT METHOD
% -----------------------------------------------------------------------------------------
% Input variables
%   NNOD    : number of nodes
%   NBC     : number of beam-column elements including truss elements
%   NT3     : number of constant strain triangular elements (T3)
%   NQ4     : number of quadrilateral elements (Q4)
%   NMAT    : numberofmaterialtypes
%   NSEC    : number of section types
%   IPLSTR  : indicator for plane stress ( = 1) or plane strain ( = 2)
%   NGAUSS  : number of 1D Gaussian points
%   NDOUBGS : number of 2D Gaussian points
%   COOR    : nodal coordinates
%   NFIX    : boundary conditions for nodes
%   EXLD    : external load
%   VECTY   : direction of the local y-axis of beam-cloumn elements
%   IDBC    : identification matrix for beam-column elements
%   IDT3    : identification matrix for T3 elements
%   IDQ4    : identification matrix for Q4 elements
%   PROP    : material properties
%   SECT    : section properties
% -----------------------------------------------------------------------------------------
%                                   PART 1 OF FEM PROGRAM                                                                  
% -----------------------------------------------------------------------------------------
clc;

clear;
filename= input('Please input the filename: ex. test.ipt\n','s');
infile= fopen(filename,'r');

cellarray= textscan(infile,'%s');

title1= cellarray{1}{1};
title2= cellarray{1}{2};
funit= cellarray{1}{3};
lunit= cellarray{1}{4};

NNOD= str2double(cellarray{1}{14});
NBC= str2double(cellarray{1}{15});
NT3= str2double(cellarray{1}{16});
NQ4= str2double(cellarray{1}{17});
NMAT= str2double(cellarray{1}{18});
NSEC= str2double(cellarray{1}{19});
IPLSTR= str2double(cellarray{1}{20});
NGAUSS= str2double(cellarray{1}{21});
NDOUBGS= NGAUSS*NGAUSS;

[COOR,NFIX,EXLD,VECTY,IDBC,IDT3,IDQ4,PROP,SECT]=inputFEM(NNOD,NBC,NT3,NQ4,NMAT,NSEC,cellarray);

% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
%                                   PART 2 OF FEM PROGRAM
% -----------------------------------------------------------------------------------------
[NEQ,IDND]= idmat(NNOD,NFIX);
[MHB,MSS,LMBC,LMT3,LMQ4,NONZ,NADD]= skyline(NBC,NT3,NQ4,NEQ,IDND,IDBC,IDT3,IDQ4);
[GLOAD]= load(NEQ,NNOD,IDND,EXLD);


% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
%                                   PART 3 OF FEM PROGRAM
% -----------------------------------------------------------------------------------------
SS= zeros(MSS,1);
[SS]= stiffbc(NBC,COOR,VECTY,IDBC,PROP,SECT,NADD,LMBC,SS);

% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
%                                   PART 4 OF FEM PROGRAM
% -----------------------------------------------------------------------------------------
[SS]= stifft3(NT3,COOR,IDT3,PROP,SECT,LMT3,NADD,IPLSTR,SS);

% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
%                                   PART 5 OF FEM PROGRAM
% -----------------------------------------------------------------------------------------
if(NQ4~= 0)
  [PLA,WGT]= gauss(NGAUSS);
  [SS,NDOUBGS]= stiffq4(NQ4,COOR,SS,IDQ4,PROP,SECT,PLA,WGT,NGAUSS,IPLSTR,NADD,LMQ4);
end
[SS,ISIGN]= decom(NONZ,NADD,SS,MHB,NEQ);
if(ISIGN == 1)
  disp('[Error]: Zero diagonal elements in SS matrix ');
  return;
end
NLC=1;
[GDISP]= solver(SS,GLOAD,NEQ,MHB,NLC,NONZ,NADD);

% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
%                                   PART 6 OF FEM PROGRAM
% -----------------------------------------------------------------------------------------
if(NBC~= 0)
  [ELFOR]= forcebc(NBC,COOR,IDBC,VECTY,PROP,SECT,LMBC,GDISP);
end
if(NT3~= 0)
  [STRAINT3,STRESST3]= forcet3(IPLSTR,NT3,COOR,IDT3,PROP,LMT3,GDISP);
end
if(NQ4~= 0)
[STRAINQ4,STRESSQ4,COORGAUS]= forceq4(NGAUSS,NDOUBGS,NQ4,IPLSTR,COOR,IDQ4,PROP,LMQ4,GDISP);
end

% -----------------------------------------------------------------------------------------
%                                     END OF FEM PROGRAM
% -----------------------------------------------------------------------------------------



