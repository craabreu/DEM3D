function [] = readfile(name)

global dim LengthScale TimeScale Time Lbox ...
       NC Component_Number Component_Radius ...
       NW Wall_IndexType Wall_L1 Wall_L2 Wall_L3 Wall_Position Wall_U Wall_W Wall_N ...
       NP Particle_TypeIndex Particle_Position Vectors Particle_U Particle_W Particle_N

global Particle_Velocity

file = fopen(name,'r');

ncc = 2;
huge = 2^(8*ncc-1)-2;
Complete = 72;
A = 1.25*pi/huge;
B = 0.75*pi;

ib = 4;
rb = 8;

ncc = ['int' num2str(8*ncc)];
ib = ['int' num2str(8*ib)];
rb = ['float' num2str(8*rb)];

Signature = fread(file,1,'int32');
if ~isequal(Signature,1129137476), disp('Error: Invalid configuration file.'); return; end
Version = fread(file,1,'int8');
dim = fread(file,1,'int8');
dav = 2*dim-3;
ini = 7 - 2*dim;
Form = fread(file,1, 'int8');

Full = Form == Complete;

LengthScale = fread(file,1,rb);
TimeScale   = fread(file,1,rb);
if Full, MassScale = fread(file,1,rb); end
Time = fread(file,1,rb);
Lbox = fread(file,dim,rb);

Factor = 0.5*max(Lbox)/huge;

NC = fread(file,1,ib);
Component_Number = zeros(NC,1);
Component_Radius = zeros(NC,1);
if Full
    Component_Density = zeros(NC,1);
    Component_Stiffness = zeros(NC,1);
    Component_ShearStiffness = zeros(NC,1);
    Component_Mass = zeros(NC,1);
    Component_Inertia = zeros(NC,1);
end
for i = 1:NC
    Component_Number(i) = fread(file,1,ib);
    Component_Radius(i) = fread(file,1,rb);
    if Full
        Component_Density(i) = fread(file,1,rb);
        Component_Stiffness(i) = fread(file,1,rb);
        if Version == 2
            Component_ShearStiffness(i) = fread(file,1,rb);
        else
            Component_ShearStiffness(i) = 2./7.*Component_Stiffness(i);
        end
        Component_Mass(i) = fread(file,1,rb);
        Component_Inertia(i) = fread(file,1,rb);
    end
end

NW = fread(file,1,ib);
if NW > 0
    Wall_TypeIndex = zeros(NW,1);
    Wall_L1 = zeros(NW,1);
    Wall_L2 = zeros(NW,1);
    Wall_L3 = zeros(NW,1);
    Wall_U  = zeros(NW,3);
    Wall_W  = zeros(NW,3);
    Wall_N  = zeros(NW,3);
    Wall_Position = zeros(NW,3);
    if Full
        Wall_Velocity = zeros(NW,3);
        Wall_Omega = zeros(NW,3);
        Wall_Stiffness = zeros(NW,3);
        Wall_ShearStiffness = zeros(NW,3);
    end
end
for i = 1:NW
    Wall_IndexType(i) = fread(file,1,ib);
    Wall_L1(i) = fread(file,1,rb);
    Wall_L2(i) = fread(file,1,rb);
    Wall_L3(i) = fread(file,1,rb);
    if Full
        Wall_U(i,1:dim) = fread(file,dim,rb)';
        Wall_W(i,1:dim) = fread(file,dim,rb)';
        Wall_N(i,1:dim) = fread(file,dim,rb)';
    else
        phi = fread(file,dim-1,ncc)';
        Wall_U(i,1:dim) = SphericToCart( A*phi + B );
        phi = fread(file,dim-1,ncc)';
        Wall_W(i,1:dim) = SphericToCart( A*phi + B );
        [Wall_U(i,:),Wall_W(i,:)] = GramShmidt( Wall_U(i,:),Wall_W(i,:) );
        Wall_N(i,:) = cross(Wall_U(i,:),Wall_W(i,:));
    end
    Wall_Position(i,1:dim) = fread(file,dim,rb)';
    if Full
        Wall_Velocity(i,1:dim) = fread(file,dim,rb)';
        Wall_Omega(i,ini:3) = fread(file,dav,rb)';
        if Version == 2 
            Wall_Stiffness(i) = fread(file,1,rb);
            Wall_ShearStiffness(i) = fread(file,1,rb);
        end
    end
end

NP = sum(Component_Number(1:NC));
Vectors = fread(file,1,'int8') ~= 0;

if NP > 0
    Particle_Position = zeros(NP,3);
    Particle_TypeIndex = zeros(NP,1);
    if Full
        Particle_Velocity = zeros(NP,3);
        Particle_Omega = zeros(NP,3);
    end
    if Vectors ~= 0
        Particle_U = zeros(NP,3);
        Particle_W = zeros(NP,3);
        Particle_N = zeros(NP,3);
    end
end

if Full
    for i = 1:NP
        Particle_TypeIndex(i)  = fread(file,1,ib);
        Particle_Position(i,1:dim) = fread(file,dim,rb)';
        Particle_Velocity(i,1:dim) = fread(file,dim,rb)';
        Particle_Omega(i,ini:3)    = fread(file,dav,rb)';
        if Vectors
            Particle_U(i,1:dim) = fread(file,dim,rb)';
            Particle_W(i,1:dim) = fread(file,dim,rb)';
            Particle_N(i,:) = fread(file,dim,rb)';
        end
    end
else
    i = 0;
    for comp = 1:NC
        Component_Number(comp);
        for j = 1:Component_Number(comp)
            i = i + 1;
            Particle_TypeIndex(i) = comp;
            Particle_Position(i,1:dim) = Factor*fread(file,dim,ncc)';
            if Vectors
                phi = fread(file,dim-1,ncc)';
                Particle_U(i,1:dim) = SphericToCart( A*phi + B );
                phi = fread(file,dim-1,ncc)';
                Particle_W(i,1:dim) = SphericToCart( A*phi + B );
                [Particle_U(i,:),Particle_W(i,:)] = GramShmidt( Particle_U(i,:),Particle_W(i,:) );
                Particle_N(i,:) = cross(Particle_U(i,:),Particle_W(i,:));
            end
        end
    end
end

file = fclose(file);

Lbox = LengthScale*Lbox;
Component_Radius = LengthScale*Component_Radius;
Wall_L1 = LengthScale*Wall_L1;
Wall_L2 = LengthScale*Wall_L2;
Wall_L3 = LengthScale*Wall_L3;
Wall_Position = LengthScale*Wall_Position;
Particle_Position = LengthScale*Particle_Position;

% Component_Number
% Component_Radius
% [Particle_TypeIndex Particle_Position]

%----------------------------------------------------------------------

function U = SphericToCart( phi )

n = size(phi);
if n == 1
    U = [cos(phi) sin(phi)];
else
    U = [ cos(phi(2))*[cos(phi(1)) sin(phi(1))] sin(phi(2))];
end

%----------------------------------------------------------------------

function [U,V] = GramShmidt( u, v )

U = u/norm(u);
V = v - dot(U,v)*U;
V = V/norm(V);

%----------------------------------------------------------------------
