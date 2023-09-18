% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: inputFEM
% -----------------------------------------------------------------------------------------
% Input variables
%   NNOD         : number of nodes
%   NBC          : number of beam-column elements including truss elements
%   NT3          : number of constant strain triangular elements (T3)
%   NQ4          : number of quadrilateral elements (Q4)
%   NMAT         : numberofmaterialtypes
%   NSEC         : number of section types
%   IPLSTR       : indicator for plane stress ( = 1) or plane strain ( = 2)
%   NGAUSS       : number of 1D Gaussian points
%   NDOUBGS      : number of 2D Gaussian points
%   COOR(3,NNOD) : nodal coordinates
%     COOR(1,*)  : x-axis coordinates
%     COOR(2,*)  : y-axis coordinates
%     COOR(3,*)  : z-axis coordinates
%   NFIX(6,NNOD)    : boundary conditions for nodes
%     NFIX(1,*)  : -1 if fixed, 0 if free for u
%     NFIX(2,*)  : -1 if fixed, 0 if free for v
%     NFIX(3,*)  : -1 if fixed, 0 if free for w
%     NFIX(4,*)  : -1 if fixed, 0 if free for theta_x
%     NFIX(5,*)  : -1 if fixed, 0 if free for theta_y
%     NFIX(6,*)  : -1 if fixed, 0 if free for theta_z
%   EXLD(6,NNOD)    : external load
%     EXLD(1,*)  : force of X-direction
%     EXLD(2,*)  : force of Y-direction
%     EXLD(3,*)  : force of Z-direction
%     EXLD(4,*)  : moment of X-direction
%     EXLD(5,*)  : moment of Y-direction
%     EXLD(6,*)  : moment of Z-direction
%   VECTY(3,NBC)   : direction of the local y-axis of beam-cloumn elements
%     VECTY(1,*) : X-direction of the local y-axis of BC elements
%     VECTY(2,*) : Y-direction of the local y-axis of BC elements
%     VECTY(3,*) : Z-direction of the local y-axis of BC elements
%   IDBC(6,NBC)    : identification matrix for beam-column elements
%     IDBC(1,*)  : node1
%     IDBC(2,*)  : node2
%     IDBC(3,*)  : node3
%     IDBC(4,*)  : node4
%     IDBC(5,*)  : blank
%     IDBC(6,*)  : blank
%   IDT3(6,NT3)    : identification matrix for T3 elements
%     IDT3(1,*)  : node1
%     IDT3(2,*)  : node2
%     IDT3(3,*)  : node3
%     IDT3(4,*)  : material
%     IDT3(5,*)  : section
%     IDT3(6,*)  : blank
%   IDQ4(6,NQ4)    : identification matrix for Q4 elements
%     IDQ4(1,*)  : node1
%     IDQ4(2,*)  : node2
%     IDQ4(3,*)  : node3
%     IDQ4(4,*)  : node4
%     IDQ4(5,*)  : material
%     IDQ4(6,*)  : section
%   PROP(4,NMAT)    : material properties
%     PROP(1,*)  : Young's modulus
%     PROP(2,*)  : Poisson's ratio
%     PROP(3,*)  : blank
%     PROP(4,*)  : blank
%   SECT(6,NSEC) : section properties
%     SECT(1,*)  : cross-sectional area A(beam-column element) or thinkness t(T3 or Q4 element)
%     SECT(2,*)  : I of y-axis
%     SECT(3,*)  : I of Z-axis
%     SECT(4,*)  : J
%     SECT(5,*)  : blank
%     SECT(6,*)  : blank
% -----------------------------------------------------------------------------------------
function[COOR,NFIX,EXLD,VECTY,IDBC,IDT3,IDQ4,PROP,SECT]=inputFEM(NNOD,NBC,NT3,NQ4,NMAT,NSEC,cellarray)

  cellno= 23;
  COOR= zeros(3,NNOD);
  for m= 1:NNOD  
    cellno= cellno+1;  
    for n= 1:3
      cellno= cellno+1;
      COOR(n,m)= str2double(cellarray{1}{cellno});
    end
  end

  cellno= cellno+2;
  NFIX= zeros(6,NNOD);
  for m= 1:NNOD  
    cellno= cellno+1;  
    for n= 1:6 
      cellno= cellno+1;
      NFIX(n,m)= str2double(cellarray{1}{cellno});
    end
  end

  cellno= cellno+2;
  EXLD= zeros(6,NNOD);
  for m= 1:NNOD  
    cellno= cellno+1;  
    for n= 1:6 
      cellno= cellno+1;
      EXLD(n,m)= str2double(cellarray{1}{cellno});
    end
  end

  cellno= cellno+2;
  VECTY= zeros(3,NBC);
  for m= 1:NBC
    cellno= cellno+1;
    for n= 1:3
      cellno= cellno+1;
      VECTY(n,m)=str2double(cellarray{1}{cellno});
    end
  end

  cellno= cellno+2;
  IDBC= zeros(6,NBC);
  for m= 1:NBC
    cellno= cellno+1;
    for n=1:6
      cellno= cellno+1;
      IDBC(n,m)=str2double(cellarray{1}{cellno});
    end
  end

  cellno= cellno+2;
  IDT3= zeros(6,NT3);
  for m= 1:NT3
    cellno= cellno+1;  
    for n= 1:6  
      cellno= cellno+1;
      IDT3(n,m)= str2double(cellarray{1}{cellno});
    end
  end

  cellno= cellno+2;
  IDQ4= zeros(6,NQ4);
  for m= 1:NQ4 
    cellno= cellno+1;  
    for n= 1:6  
      cellno= cellno+1;
      IDQ4(n,m)= str2double(cellarray{1}{cellno});
    end
  end


  cellno= cellno+2;
  PROP= zeros(4,NMAT);
  for m= 1:NMAT 
    cellno= cellno+1;  
    for n= 1:4  
      cellno= cellno+1;
      PROP(n,m)= str2double(cellarray{1}{cellno});
    end
  end

  cellno= cellno+2;
  SECT= zeros(6,NSEC);
  for m= 1:NSEC 
    cellno= cellno+1;  
    for n= 1:6  
      cellno= cellno+1;
      SECT(n,m)= str2double(cellarray{1}{cellno});
    end
  end
end
% -----------------------------------------------------------------------------------------