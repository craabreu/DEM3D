function [] = PlotConfig( filename, WallColor, Color, PR, WR, Opacity, MyView, Half, PBC, MaxNW )

global dim LengthScale TimeScale Time Lbox ...
       NC Component_Number Component_Radius ...
       NW Wall_IndexType Wall_L1 Wall_L2 Wall_L3 Wall_Position Wall_U Wall_W Wall_N ...
       NP Particle_TypeIndex Particle_Position Vectors Particle_U Particle_W Particle_N

if nargin < 6
    disp('The function PlotConfig requires at least 6 arguments');
    return;
end
if nargin < 10
    MaxNW = NW;
    if nargin < 9
        PBC = [0 0 0];
        if nargin < 8
            Half = 0;
            if nargin < 7
                MyView = 3;
            end
        end
    end
end

% Read the configuration from the file:
readfile(filename);

% Draw walls:
if MaxNW < 0, MaxNW = NW; end
for i = 1:MaxNW
    type = Wall_IndexType(i);
    if type == 1
        [x y z] = rect(2*Wall_L1(i), 2*Wall_L2(i), Wall_Position(i,:), Wall_U(i,:), Wall_W(i,:));
    elseif type == 2
        [x y z] = cone(WR, 2*Wall_L1(i), Wall_L2(i), Wall_L3(i), Wall_Position(i,:), Wall_U(i,:), Wall_W(i,:), Wall_N(i,:), Half );
    end
    p = surf2patch(x,y,z);
    patch(p,'facecolor',WallColor,'edgecolor','none');
end

% Set the wall opacity:
alpha(Opacity);

if Vectors
    % Draw the particles:
    for i = 1:NP
        if ~Half | (Half & Particle_Position(i,1) > 0)
            comp = Particle_TypeIndex(i);
            [x,y,z] = sphere( PR, Component_Radius(comp), Particle_Position(i,:), ...
                Particle_U(i,:), Particle_W(i,:), Particle_N(i,:) );
            p = surf2patch(x,y,z);
            patch(p,'facecolor',Color(comp,:),'edgecolor','none');
        end
    end
else
    % Construct a model for the particles:
    [xe,ye,ze] = sphere(PR,1,[0 0 0],[1 0 0],[0 1 0],[0 0 1]);
    % Draw the particles:
    for i = 1:NP
        if ~Half | (Half & Particle_Position(i,1) > 0)
            r    = Particle_Position(i,:);
            comp = Particle_TypeIndex(i);
            Rad  = Component_Radius(comp);
            p    = surf2patch(Rad*xe+r(1),Rad*ye+r(2),Rad*ze+r(3));
            patch(p,'facecolor',Color(comp,:),'edgecolor','none');
        end
    end
end

% Set figure properties:
L = Lbox + 2.*PBC'*max(Component_Radius);
axis(0.5*[-L(1) L(1) -L(2) L(2) -L(3) L(3)]);
axis off;
daspect([1 1 1]);
view(MyView);
lighting gouraud;
camlight;
