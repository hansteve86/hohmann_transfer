%% RK4 HOHMANN Animation MATLAB SCRIPT 
clc;clear;close all

% Parameters
Me = 1.898e27;% mass of body a
G = 6.67408e-11;    % gravitational constant
Mu = G*Me;
h = 50; %time step (seconds)

%equations of motion written in state space form
f1 = @(t, r, rdot, theta, thetad) rdot; 
f2 = @(t, r, rdot, theta, thetad) -Mu/r^2 + r*thetad^2;
f3 = @(t, r, rdot, theta, thetad) thetad;
f4 = @(t, r, rdot, theta, thetad) -2*rdot*thetad/r;

% Initial Conditions (testing random ICs)

r(1) = 71492000; %inital position of satellite distance from planet 1 to planet 2)
r2(1) = 91492000; %inital position of satellite
a = (r(1) + r2(1))/2;
rdot(1) = 0; % circular orbit
theta(1) = pi+0.5*pi; % circular orbit
thetad(1) = sqrt(Mu/r(1)^3); %circular orbit
tf = sqrt(4*pi^2*r(1)^3/Mu);

% initial time
t(1) = 0;
%% inner orbit
for i = 1:tf/h
    
    t(i+1) = t(i) + h;
    
    k1r = f1(t(i), r(i), rdot(i), theta(i), thetad(i));
    k1rdot = f2(t(i), r(i), rdot(i), theta(i), thetad(i));
    k1theta = f3(t(i), r(i), rdot(i), theta(i), thetad(i));
    k1thetad = f4(t(i), r(i), rdot(i), theta(i), thetad(i));
     
    k2r = f1(t(i)+h/2, r(i)+k1r*h/2, rdot(i)+k1rdot*h/2, theta(i)+k1theta*h/2, thetad(i)+k1thetad*h/2);
    k2rdot = f2(t(i)+h/2, r(i)+k1r*h/2, rdot(i)+k1rdot*h/2, theta(i)+k1theta*h/2, thetad(i)+k1thetad*h/2);
    k2theta = f3(t(i)+h/2, r(i)+k1r*h/2, rdot(i)+k1rdot*h/2, theta(i)+k1theta*h/2, thetad(i)+k1thetad*h/2);
    k2thetad = f4(t(i)+h/2, r(i)+k1r*h/2, rdot(i)+k1rdot*h/2, theta(i)+k1theta*h/2, thetad(i)+k1thetad*h/2);
    
    k3r = f1(t(i)+h/2, r(i)+k2r*h/2, rdot(i)+k2rdot*h/2, theta(i)+k2theta*h/2, thetad(i)+k2thetad*h/2);
    k3rdot = f2(t(i)+h/2, r(i)+k2r*h/2, rdot(i)+k2rdot*h/2, theta(i)+k2theta*h/2, thetad(i)+k2thetad*h/2);
    k3theta = f3(t(i)+h/2, r(i)+k2r*h/2, rdot(i)+k2rdot*h/2, theta(i)+k2theta*h/2, thetad(i)+k2thetad*h/2);
    k3thetad = f4(t(i)+h/2, r(i)+k2r*h/2, rdot(i)+k2rdot*h/2, theta(i)+k2theta*h/2, thetad(i)+k2thetad*h/2);
   
    
    k4r = f1(t(i)+h, r(i)+k3r*h, rdot(i)+k3rdot*h, theta(i)+k3theta*h, thetad(i)+k3thetad*h);
    k4rdot = f2(t(i)+h, r(i)+k3r*h, rdot(i)+k3rdot*h, theta(i)+k3theta*h, thetad(i)+k3thetad*h);
    k4theta = f3(t(i)+h, r(i)+k3r*h, rdot(i)+k3rdot*h, theta(i)+k3theta*h, thetad(i)+k3thetad*h);
    k4thetad = f4(t(i)+h, r(i)+k3r*h, rdot(i)+k3rdot*h, theta(i)+k3theta*h, thetad(i)+k3thetad*h);
    
    %updating position and speed
    r(i+1,1) = r(i,1) + (h/6)*(k1r+2*k2r+2*k3r+k4r);
    rdot(i+1,1) = rdot(i,1) + (h/6)*(k1rdot+2*k2rdot+2*k3rdot+k4rdot);
    theta(i+1,1) = theta(i,1) + (h/6)*(k1theta+2*k2theta+2*k3theta+k4theta);
    thetad(i+1,1) = thetad(i,1) + (h/6)*(k1thetad+2*k2thetad+2*k3thetad+k4thetad);
end




%% transfer orbit
% Parameters
tft = pi*sqrt(((r(1)+r2(1))/2)^3/(Mu)); % time take for orbit

%equations of motion written in state space form
ft1 = @(tt, rt, rdott, thetat, thetadt) rdott; 
ft2 = @(tt, rt, rdott, thetat, thetadt) -Mu/rt^2 + rt*thetadt^2; 
ft3 = @(tt, rt, rdott, thetat, thetadt) thetadt; 
ft4 = @(tt, rt, rdott, thetat, thetadt) -2*rdott*thetadt/rt; 

% Initial Conditions

rt(1) = r(1); %inital position of satellite, moment before first burn
NumRDotEntries = length(rdot);
rdott(1) = rdot(NumRDotEntries); % transfer ellipitcal orbit, final v value from inner orbit + delta v
thetat(1) = theta(1); % circular orbit
thetadt(1) = sqrt(Mu/r(NumRDotEntries)^3) + ((sqrt(Mu/r(1))*(sqrt(2*r2(1)/(r(1)+r2(1)))-1))/r(1)) +  ((sqrt(Mu/r2(1))*(sqrt(2*r(1)/(r(1)+r2(1)))+1))/r2(1))/r2(1); %circular orbit angular velocity

% initial time)
tt(1) = tft+max(t); %time for orbit + end time of initial orbit

%% transfer orbit rk4
for i = 1:tft/h
    
    tt(i+1) = tt(i) + h;
    
    k1rt = ft1(tt(i), rt(i), rdott(i), thetat(i), thetadt(i));
    k1rdott = ft2(tt(i), rt(i), rdott(i), thetat(i), thetadt(i));
    k1thetat = ft3(tt(i), rt(i), rdott(i), thetat(i), thetadt(i));
    k1thetadt = ft4(tt(i), rt(i), rdott(i), thetat(i), thetadt(i));
     
    k2rt = ft1(tt(i)+h/2, rt(i)+k1rt*h/2, rdott(i)+k1rdott*h/2, thetat(i)+k1thetat*h/2, thetadt(i)+k1thetadt*h/2);
    k2rdott = ft2(tt(i)+h/2, rt(i)+k1rt*h/2, rdott(i)+k1rdott*h/2, thetat(i)+k1thetat*h/2, thetadt(i)+k1thetadt*h/2);
    k2thetat = ft3(tt(i)+h/2, rt(i)+k1rt*h/2, rdott(i)+k1rdott*h/2, thetat(i)+k1thetat*h/2, thetadt(i)+k1thetadt*h/2);
    k2thetadt = ft4(tt(i)+h/2, rt(i)+k1rt*h/2, rdott(i)+k1rdott*h/2, thetat(i)+k1thetat*h/2, thetadt(i)+k1thetadt*h/2);
    
    k3rt = ft1(tt(i)+h/2, rt(i)+k2rt*h/2, rdott(i)+k2rdott*h/2, thetat(i)+k2thetat*h/2, thetadt(i)+k2thetadt*h/2);
    k3rdott = ft2(tt(i)+h/2, rt(i)+k2rt*h/2, rdott(i)+k2rdott*h/2, thetat(i)+k2thetat*h/2, thetadt(i)+k2thetadt*h/2);
    k3thetat = ft3(tt(i)+h/2, rt(i)+k2rt*h/2, rdott(i)+k2rdott*h/2, thetat(i)+k2thetat*h/2, thetadt(i)+k2thetadt*h/2);
    k3thetadt = ft4(tt(i)+h/2, rt(i)+k2rt*h/2, rdott(i)+k2rdott*h/2, thetat(i)+k2thetat*h/2, thetadt(i)+k2thetadt*h/2);
   
    
    k4rt = ft1(tt(i)+h, rt(i)+k3rt*h, rdott(i)+k3rdott*h, thetat(i)+k3thetat*h, thetadt(i)+k3thetadt*h);
    k4rdott = ft2(tt(i)+h, rt(i)+k3rt*h, rdott(i)+k3rdott*h, thetat(i)+k3thetat*h, thetadt(i)+k3thetadt*h);
    k4thetat = ft3(tt(i)+h, rt(i)+k3rt*h, rdott(i)+k3rdott*h, thetat(i)+k3thetat*h, thetadt(i)+k3thetadt*h);
    k4thetadt = ft4(tt(i)+h, rt(i)+k3rt*h, rdott(i)+k3rdott*h, thetat(i)+k3thetat*h, thetadt(i)+k3thetadt*h);
    
    %updating position and speed
    rt(i+1,1) = rt(i,1) + (h/6)*(k1rt+2*k2rt+2*k3rt+k4rt);
    rdott(i+1,1) = rdott(i,1) + (h/6)*(k1rdott+2*k2rdott+2*k3rdott+k4rdott);
    thetat(i+1,1) = thetat(i,1) + (h/6)*(k1thetat+2*k2thetat+2*k3thetat+k4thetat);
    thetadt(i+1,1) = thetadt(i,1) + (h/6)*(k1thetadt+2*k2thetadt+2*k3thetadt+k4thetadt);
end





%% outer orbit

f5 = @(t, r2, rdot2, theta2, thetad2) rdot2;
f6 = @(t, r2, rdot2, theta2, thetad2) -Mu/r2^2 + r2*thetad2^2;
f7 = @(t, r2, rdot2, theta2, thetad2) thetad2;
f8 = @(t, r2, rdot2, theta2, thetad2) -2*rdot2*thetad2/r2;

rdot2(1) = 0; % circular orbit
theta2(1) = pi/2; % circular orbit
thetad2(1) = sqrt(Mu/r2(1)^3); %circular orbit
t2(1) = tf;
tf2 = sqrt(4*pi^2*r2(1)^3/Mu); %start from the end time of transfer orbit (add the duration of the orbit onto this)
h2 = h;
for i = 1:tf2/h2
    
    t2(i+1) = t2(i) + h2;
    
    k1r2 = f5(t2(i), r2(i), rdot2(i), theta2(i), thetad2(i));
    k1rdot2 = f6(t2(i), r2(i), rdot2(i), theta2(i), thetad2(i));
    k1theta2 = f7(t2(i), r2(i), rdot2(i), theta2(i), thetad2(i));
    k1thetad2 = f8(t2(i), r2(i), rdot2(i), theta2(i), thetad2(i));
     
    k2r2 = f5(t2(i)+h2/2, r2(i)+k1r2*h2/2, rdot2(i)+k1rdot2*h2/2, theta2(i)+k1theta2*h2/2, thetad2(i)+k1thetad2*h2/2);
    k2rdot2 = f6(t2(i)+h2/2, r2(i)+k1r2*h2/2, rdot2(i)+k1rdot2*h2/2, theta2(i)+k1theta2*h2/2, thetad2(i)+k1thetad2*h2/2);
    k2theta2 = f7(t2(i)+h2/2, r2(i)+k1r2*h2/2, rdot2(i)+k1rdot2*h2/2, theta2(i)+k1theta2*h2/2, thetad2(i)+k1thetad2*h2/2);
    k2thetad2 = f8(t2(i)+h2/2, r2(i)+k1r2*h2/2, rdot2(i)+k1rdot2*h2/2, theta2(i)+k1theta2*h2/2, thetad2(i)+k1thetad2*h2/2);
    
    k3r2 = f5(t2(i)+h2/2, r2(i)+k2r2*h2/2, rdot2(i)+k2rdot2*h2/2, theta2(i)+k2theta2*h2/2, thetad2(i)+k2thetad2*h2/2);
    k3rdot2 = f6(t2(i)+h2/2, r2(i)+k2r2*h2/2, rdot2(i)+k2rdot2*h2/2, theta2(i)+k2theta2*h2/2, thetad2(i)+k2thetad2*h2/2);
    k3theta2 = f7(t2(i)+h2/2, r2(i)+k2r2*h2/2, rdot2(i)+k2rdot2*h2/2, theta2(i)+k2theta2*h2/2, thetad2(i)+k2thetad2*h2/2);
    k3thetad2 = f8(t2(i)+h2/2, r2(i)+k2r2*h2/2, rdot2(i)+k2rdot2*h2/2, theta2(i)+k2theta2*h2/2, thetad2(i)+k2thetad2*h2/2);
   
    
    k4r2 = f5(t2(i)+h2, r2(i)+k3r2*h2, rdot2(i)+k3rdot2*h2, theta2(i)+k3theta2*h2, thetad2(i)+k3thetad2*h2);
    k4rdot2 = f6(t2(i)+h2, r2(i)+k3r2*h2, rdot2(i)+k3rdot2*h2, theta2(i)+k3theta2*h2, thetad2(i)+k3thetad2*h2);
    k4theta2 = f7(t2(i)+h2, r2(i)+k3r2*h2, rdot2(i)+k3rdot2*h2, theta2(i)+k3theta2*h2, thetad2(i)+k3thetad2*h2);
    k4thetad2 = f8(t2(i)+h2, r2(i)+k3r2*h2, rdot2(i)+k3rdot2*h2, theta2(i)+k3theta2*h2, thetad2(i)+k3thetad2*h2);
    
    %updating position and speed
    r2(i+1,1) = r2(i,1) + (h2/6)*(k1r2+2*k2r2+2*k3r2+k4r2);
    rdot2(i+1,1) = rdot2(i,1) + (h2/6)*(k1rdot2+2*k2rdot2+2*k3rdot2+k4rdot2);
    theta2(i+1,1) = theta2(i,1) + (h2/6)*(k1theta2+2*k2theta2+2*k3theta2+k4theta2);
    thetad2(i+1,1) = thetad2(i,1) + (h2/6)*(k1thetad2+2*k2thetad2+2*k3thetad2+k4thetad2);
end

%% Animation
FPS = 10;

Video = VideoWriter("Animation1");
%Video.Quality = 100;
Video.FrameRate = FPS;
open(Video);  

ShipScale = 0.1;
MainBodyScale = 0.2;
SatelliteScale = 0.07;

FigureMaxX = max(r2.*cos(theta2));
FigureMaxY = max(r2.*cos(theta2));
FigureScale = max([FigureMaxX FigureMaxY]);

PrevShipX = 0;
PrevShipY = 0;

figure(1)

InnerRadii = r;
InnerThetas = theta;

HohmannRadii = rt;
HohmannThetas = thetat;

OuterRadii = r2;
OuterThetas = theta2;

for n = 1:(length(t) + length(tt) + length(t2))
    %% Setup plot
    clf; %clear the plot
    axis([-FigureMaxX FigureMaxX -FigureMaxY FigureMaxY]*1.15) %expand our axes by 15%
    %point format is [x1,x2],[y1,y2]
    
    Index = n;
    
    if n > length(t) + length(tt) % change theta according to our progression
        Index = Index - length(t) - length(tt);
        Radius = OuterRadii(Index); 
        Theta = OuterThetas(Index); 
    elseif n > length(t)
        Index = Index - length(t);
        Radius = HohmannRadii(Index);
        Theta = HohmannThetas(Index);   
    else
        Radius = InnerRadii(Index);
        Theta = InnerThetas(Index);            
    end
    
    xlabel('x-coordinates (m)');
    ylabel('y-coordinates (m)');
%   title(TitleString);

    hold on
    [MainBodyImgRef, ~, MainBodyAlpha] = imread('Jupiter.png');
    [SatelliteImgRef, ~, SatelliteAlpha] = imread('Earth.png');
    [ShipImgRef, ~, ShipAlpha] = imread('Shuttle.png');

    hold on
    InnerX = InnerRadii.*cos(InnerThetas);
    InnerY = InnerRadii.*sin(InnerThetas);
    plot(InnerX,InnerY,'-b','LineWidth',2); % MainBody's orbit

    hold on
    HohmannX = HohmannRadii.*cos(HohmannThetas);
    HohmannY = HohmannRadii.*sin(HohmannThetas);
    plot(HohmannX,HohmannY,'g--','LineWidth',2); %Hohmann orbit

    hold on
    OuterX = OuterRadii.*cos(OuterThetas);
    OuterY = OuterRadii.*sin(OuterThetas);
    plot(OuterX,OuterY,'-r','LineWidth',2); %Outer orbit

    % Scalar values used for proper image scaling (we chose to define image
    % sizes with a value between 0-1 (1 being as being as the figure)
    XLim = get(gca,'XLim'); 
    YLim = get(gca,'YLim');

    hold on % MainBody image
    image(MainBodyImgRef,'AlphaData', MainBodyAlpha,'XData',[-1 1]*FigureScale*MainBodyScale,'YData',[1 -1]*FigureScale*MainBodyScale);

%     hold on % Satellite image
%     SatelliteX = FigureScale * -0.4; %TEMP
%     SatelliteY = FigureScale * 0.4; %TEMP
%     image(SatelliteImgRef,'AlphaData', MainBodyAlpha,'XData',SatelliteX + [-1 1]*FigureScale*SatelliteScale,'YData',SatelliteY + [1 -1]*FigureScale*SatelliteScale);

    hold on % Ship image
    
    ShipX = Radius*cos(Theta);
    ShipY = Radius*sin(Theta);

    ShipHeading = rad2deg(atan2(ShipY - PrevShipY,ShipX - PrevShipX)) - 90; %get the heading locally off of the previous position

    if n > 1 % Skip the first frame so we can collect an initial previous position for our heading (prevents weird 0 degree snap)
        RotatedShip = imrotate(ShipImgRef,ShipHeading,'nearest','crop');
        ShipMask = zeros(size(RotatedShip, 1), size(RotatedShip, 2), 'uint8'); % Used to make black rotated background transparent
        % RotatedShip(RotatedShip == 0) = 255; % Convert black pixels to white
        RotatedShipImg = image(RotatedShip,'AlphaData',ShipMask,'XData',ShipX + [-1 1]*FigureScale*ShipScale,'YData',ShipY + [1 -1]*FigureScale*ShipScale);
        RemainingAlpha = max(RotatedShip,[],3);% You need this to recover any non 0 values
        RotatedShipImg.AlphaData = RemainingAlpha*5; % Multiply by some constant to make the image less "ghostly"
    end
    
    PrevShipX = ShipX;
    PrevShipY = ShipY;
    
    %annotation('textbox', [0.635, 0.1, 0.1, 0.1], 'String', sprintf("%0.3f Hours",(n*h)/3600),"Color","white","FontSize",12,"LineStyle","none","LineWidth",1)
    
    hold off
    grid on
    %axis equal %this line causes the plot to shift when the shuttle moves
    axis square
    
    ax = gca; % Reference to the axes
    ax.Color = 'k';
    ax.XColor = 'r'; % Red
    ax.YColor = 'b'; % Blue
    ax.GridAlpha = 0.9;  % Grid transparency
    ax.GridColor = [0.05, 0.4, 0.1]; 
    
    CurrentFigure = gcf; %reference to the current figure
    
    Frame = getframe(CurrentFigure);%"capture" the current figure as a frame
    writeVideo(Video,Frame); %write the frame to our video
end


