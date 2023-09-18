% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: skyline
% -----------------------------------------------------------------------------------------
% Input variables
%   NBC             : number of beam-column elements including truss elements
%   NT3             : number of constant strain triangular elements 
%   NQ4             : number of quadrilateral elements
%   NEQ             : number of equations
%   IDND            : DOFs numbering
%   IDBC            : identification matrix for beam-column elements
%   IDT3            : identification matrix for T3 elements
%   IDQ4            : identification matrix for Q4 elements
% Output variables
%   MHB             : mxaimum half bandwidth (semi-band width) = max(Bi)
%   MSS             : maximum storage size for [K] in one-dimensional array
%   LMBC            : element location matrices of beam-column elements including truss elements
%   LMT3            : element location matrices of constant strain triangular elements 
%   LMQ4            : element location matrices of quadrilateral elements
%   NONZ            : first nonzero entry of each row of [K](Ii)
%   NADD            : address counter of each row of [K](Ji)
% -----------------------------------------------------------------------------------------
function[MHB,MSS,LMBC,LMT3,LMQ4,NONZ,NADD]= skyline(NBC,NT3,NQ4,NEQ,IDND,IDBC,IDT3,IDQ4)

%
% Find LMBC, LMT3, LMQ4 matrix
%
LMBC= zeros(12,NBC);
LMT3= zeros(6,NT3);
LMQ4= zeros(8,NQ4);

for m= 1:NBC
    LMBC(1:6,m)= IDND(1:6,IDBC(1,m));
    LMBC(7:12,m)= IDND(1:6,IDBC(2,m)); 
end

for m= 1:NT3
    LMT3(1:2,m)= IDND(1:2,IDT3(1,m));
    LMT3(3:4,m)= IDND(1:2,IDT3(2,m));
    LMT3(5:6,m)= IDND(1:2,IDT3(3,m));
end

for m= 1:NQ4
    LMQ4(1:2,m)= IDND(1:2,IDQ4(1,m));
    LMQ4(3:4,m)= IDND(1:2,IDQ4(2,m));
    LMQ4(5:6,m)= IDND(1:2,IDQ4(3,m));
    LMQ4(7:8,m)= IDND(1:2,IDQ4(4,m));

end

%
% Make VIBJS matrix with NEQ rows and (NBC+NT3+NQ4) columns (Chapter5, page5)
%
VIBJS= zeros(NEQ,NBC+NT3+NQ4);
for m= 1: NBC
  MINI=NEQ;
  for i=1:12
    if( LMBC(i,m)~=0 && LMBC(i,m)<MINI )
	  MINI=LMBC(i,m);
    end
  end
  for i=1:12
    if( LMBC(i,m)~=0)
      VIBJS(LMBC(i,m),m)=MINI;
    end
  end
  for i=1:NEQ
  c=VIBJS(i,:);
  t=min(c(c~=0));
    if(VIBJS(i,m)~=0 && VIBJS(i,m)>t)
      VIBJS(i,m)=t;
    end
  end
end
for m= NBC+1:NBC+NT3
  MINI=NEQ;
  for i=1:6
    if( LMT3(i,m-NBC)~=0 && LMT3(i,m-NBC)<MINI )
	  MINI=LMT3(i,m-NBC);
    end
  end
  for i=1:6
    if( LMT3(i,m-NBC)~=0)
      VIBJS(LMT3(i,m-NBC),m)=MINI;
    end
  end
  for i=1:NEQ
  c=VIBJS(i,:);
  t=min(c(c~=0));
    if(VIBJS(i,m)~=0 && VIBJS(i,m)>t)
      VIBJS(i,m)=t;
    end
  end
end
for m= NBC+NT3+1:NBC+NT3+NQ4
  MINI=NEQ;
  for i=1:8
    if( LMQ4(i,m-NBC-NT3)~=0 && LMQ4(i,m-NBC-NT3)<MINI )
	  MINI=LMQ4(i,m-NBC-NT3);
    end
  end
  for i=1:8
    if( LMQ4(i,m-NBC-NT3)~=0)
      VIBJS(LMQ4(i,m-NBC-NT3),m)=MINI;
    end
  end
  for i=1:NEQ
  c=VIBJS(i,:);
  t=min(c(c~=0));
    if(VIBJS(i,m)~=0 && VIBJS(i,m)>t)
      VIBJS(i,m)=t;
    end
  end
end

%
% Find NONZ vector (which is also named as Ii in lecture note)
%
NONZ= zeros(1,NEQ);
for m=1:NEQ
    c=VIBJS(m,:);
    NONZ(1,m)=min(c(c~=0));
end

%
% Find Bi vector
%
Bi= zeros(1,NEQ);
for m=1:NEQ
    Bi(1,m)=m-NONZ(1,m)+1;
end

% %
% % Find Si vector
% %
Si= zeros(1,NEQ);
for m=1:NEQ
  Si(1,m)=sum(Bi(1:m));
end

% %
% % Find NADD vector (which is also named as Ji in lecture note)
% %
NADD= zeros(1,NEQ);
for m=1:NEQ
  NADD(1,m)=Si(1,m)-m;
end

% %
% % Find MHB and MSS (Chapter5, page6)
% %
MHB=1;
for i=1:NEQ
    if(i-NONZ(1,i)+1 > MHB)
        MHB=i-NONZ(1,i)+1;
    end
end
MSS=NADD(NEQ)+NEQ;

end
