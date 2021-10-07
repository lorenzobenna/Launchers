function [s]=AST_pelvec(e,xmu)
%AST_PELVEC Convert from orbital elements to cartesian state vector.
%   S = AST_pelvec(E,XMU) returns the state array nx6 given the classical
%   orbital elements E(nx6) and the gravitational constant of the central
%   body XMU.  
%
%   S will be a nx6 array whose first three components, colum-wise, will be
%   the position X, Y and Z in km and the last three components will be the
%   velocity VX , VY and VZ in km/s. 
%
%   The orbital state vector must be given as:
%        E(n,1)  a    Semimajor axis    [km ]   for e ~= 1*
%                p    Orbital parameter [km ]   for e == 1
%        E(n,2)  e    Eccentricity      [-- ]       >0      
%        E(n,3)  i    Inclination       [rad]   0 < i < pi  
%        E(n,4)  o    Ascending node    [rad]   0 < o < 2*pi
%        E(n,5)  w    Perigee Argument  [rad]   0 < w < 2*pi
%        E(n,6)  u    True Anomaly      [rad]   0 < u < 2*pi
%
%   XMU must be given in km^3/s^2.
%
%   *If the orbit is hyperbolic, the energy of the orbit will be higher than
%   0, which implies that a < 0, however sometimes is common to give the
%   absolute value of the semimajor axis for hyperbolic orbits. Therefore,
%   this function neglects the sign in the semimajor axis to simplicity.
%   
%   This function admits also a set of orbital elements as a matrix,
%   the size of the matrix must be then nx6, and each orbital element set
%   is consider to be given a row-wise. S will be then a nx6 matrix whose
%   rows will be X, Y, Z, VX, VY and VZ respectively.
%
%   See also AST_pvecle, AST_kepequ, AST_paranc
%
%   Copyright <a href="matlab: web('http://www.gmv.es','-browser')">GMV S.A.</a>, 2006

%function [s]=AST_pelvec(e,xmu)
%****************************************************************************************
% pelvec V 1.0         *             FUNCTION :  AST_pelvec                                 
%***********************             AUTHOR   :  Daniel Garcia Yarnoz (ddgy)
% COPYRIGHT GMV S.A.   *             PROJECT  :  MAST4                        
% MASE DIVISION        *             MATLAB   :  V 6.1.0 (R 12.1) 
%****************************************************************************************
%
%    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%    X                                                                         X 
%    X           PLANET ORBITAL ELEMENTS TO STATE VECTOR                       X
%    X           -             --                ---                           X
%    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%
%****************************************************************************************
%
% DESCRIPTION: Transcription of orblib's pelvec.f
%              from orbital elements to state vector (x(i), v(i)),
%              with separate input of central potential 'xmu'.
%              Taking advantadge of the matlab vectorial computation the
%              vectorial form has been added.
%
% USE        : S = PELVEC(E,XMU)
%
% HISTORY    : V 1.0: 14/06/05 -ddgy- creation (recoded from FORTRAN to MATLAB)
%              V 1.1: 04/09/06 -racg- For parabolic case e(1) is taken as p
%                                     Vectorial form added.
%
%*****************************************************************************************
%
% INPUTS     : e    (nx6) [km,-,rad,rad,rad,rad] Orbital elements:
%                    a = e(n,1)                 semi-major axis (elliptic,
%                                               and hyperbolic (a< || > 0))
%                                               with e > 1
%                    p = e(n,1)                 Orbital parameter for
%                                               parabolic orbit
%                    e = e(n,2)                 eccentricity
%                    i = e(n,3)                 inclination, in interval 0 to pi
%                    o = e(n,4)                 ascending node
%                    w = e(n,5)                 arg. of pericenter
%                    v = e(n,6)                 true anomaly
%              xmu  (1x1)   [km^3/s^2]          gravity potential of central body.  
%
% OUTPUTS    : s    (1x6)  [km,km/s]          state vector (position and velocity)
% GLOBAL     :  N/A
%
%****************************************************************************************
%
% EXTERNAL FUN :  N/A
%
% INTERNAL FUN :  N/A
%
%****************************************************************************************
%
% INPUT FILES  :  N/A
%
% OUTPUT FILES :  N/A 
%
%****************************************************************************************
      % racg, for parabolic orbits e(1) is the orbital parameter
      
      [n,m] = size( e );
      
      if m == 1
          e = e';
      end
      
      if any(abs(e(:,2)-1)) < eps
          p = e(:,1);
      else
          p = abs(e(:,1).*(1.d0 - e(:,2).^2));
      end
%     safety measure for the square root
      p = max(p,1.d-30);
      f =sqrt(abs(xmu./p));
      cv = cos(e(:,6));
      ecv = 1.d0 + e(:,2).*cv;
      
      % Distance from the focus
      r = p./ecv;
      
      u = e(:,5) + e(:,6);
      cu = cos(u);
      su = sin(u);
      co = cos(e(:,4));
      so = sin(e(:,4));
      ci = cos(e(:,3));
      si = sin(e(:,3));
      cocu = co.*cu;
      sosu = so.*su;
      socu = so.*cu;
      cosu = co.*su;
      fx = cocu - sosu.*ci;
      fy = socu + cosu.*ci;
      fz = su.*si;
      % Radial velocity
      vr = f.*e(:,2).*sin(e(:,6));
      % Tangencial velocity
      vu = f.*ecv;
      s  = [r.*fx, ...
            r.*fy, ...
            r.*fz,...
            vr.*fx-vu.*(cosu+socu.*ci),...
            vr.*fy-vu.*(sosu-cocu.*ci),...
            vr.*fz+vu.*cu.*si];

    if m == 1
        s = s';
    end

return %-----------------------------------------------------------------[MAIN FUNCTION]%