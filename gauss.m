% -----------------------------------------------------------------------------------------
%                             SUBROUTINE: gauss
% -----------------------------------------------------------------------------------------
% Input variables
%   NGAUSS       : number of GAUSS point
%
% 
% Output variables
%   PLA          : Gaussian Integration point position
%   WGT          : the correponging integration weight
% -----------------------------------------------------------------------------------------
function [PLA,WGT]= gauss(NGAUSS)
if(NGAUSS== 2)
    PLA(1)= -0.5773502692;
    PLA(2)= -PLA(1);
    WGT(1)= 1;
    WGT(2)= 1;
end

if(NGAUSS== 3)
    PLA(1)= -0.7745966692;
    PLA(2)= 0.0000000000;
    PLA(3)= -PLA(1);
    WGT(1)= 0.5555555556;
    WGT(2)= 0.8888888889;
    WGT(3)= WGT(1);
end

if(NGAUSS== 4)
    PLA(1)= -0.8611363116;
    PLA(2)= -0.3399810436;
    PLA(3)= -PLA(2);
    PLA(4)= -PLA(1);
    WGT(1)= 0.3478548451;
    WGT(2)= 0.6521451549;
    WGT(3)= WGT(2);
    WGT(4)= WGT(1);
end
end
