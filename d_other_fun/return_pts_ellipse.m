function [NX,NY]=return_pts_ellipse(Position,Axes,Eccen,PA,Color,Width,CosDec,FaceColor);
% 
% [NX,NY]=plot_ellipse(Position,Axes,Eccen,PA,Color,Width,CosDec,FaceColor);
% 
% Input  : - Position [X,Y] within current axis.
%          - [Major axis, Minor axis]
%          - Eccentricity. If given (e.g., not empty []) then Minor axis
%            is calculated from major axis and eccentricity.
%          - PA [rad], easured from X axis counterclocwise .
%          - Color, default is [0 0 0] (e.g., 'k', black).
%          - Line Width, default is 1.
%          - if cos(Dec) is given then stretch X axis by 1/cos(Dec),
%            default is 1.
%          - Ellipse face color, default is 'none'.
% Output : -  X,Y points defining the ellipse 
% Modified by DRB (NCSU) 2009 
% Original Code by: Eran O. Ofek -  March 2004
%
%--------------------------------------------------------------------------
RAD = 180./pi;

DefColor     = [0 0 0];
DefWidth     = 1;
DefCosDec    = 1;
DefFaceColor = 'none';

if (nargin==4),
   Color     = DefColor;
   Width     = DefWidth;
   CosDec    = DefCosDec;
   FaceColor = DefFaceColor;
elseif (nargin==5),
   Width     = DefWidth;
   CosDec    = DefCosDec;
   FaceColor = DefFaceColor;
elseif (nargin==6),
   CosDec    = DefCosDec;
   FaceColor = DefFaceColor;
elseif (nargin==7),
   FaceColor = DefFaceColor;
elseif (nargin==8),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(Eccen)==1),
   MajorAxis = Axes(1);
   MinorAxis = Axes(2);
else
   MajorAxis = Axes(1);
   MinorAxis = sqrt(MajorAxis(1).^2.*(1 - Eccen.^2));
end


Theta = [0:5:360]'./RAD;

X     = MajorAxis.*cos(Theta);
Y     = MinorAxis.*sin(Theta);

NX    = Position(1) + (cos(PA).*X - sin(PA).*Y)./CosDec;
NY    = Position(2) + sin(PA).*X + cos(PA).*Y;


%H     = plot(NX,NY);
H      = patch(NX,NY,'k');
set(H,'EdgeColor',Color,'FaceColor',FaceColor,'LineWidth',Width);

