function [] = Configuration( filename, PR, WR, Opacity, Skip, FPS, WallColor, Color, Width, Height, View, Half, PBC, Video, MaxNW )

if nargin == 0

    % Open a dialog box for the user to choose a file:
    [filename path] = uigetfile({'*.cfg','Configuration files (*.cfg)'},'Choose a configuration file...');
    if isequal(filename,0) | isequal(path,0)
        return;
    else
        cd(path);
        filename = [path filename];
    end

end

BaseName = GetBaseName(filename);
NDigits = GetNDigits(BaseName);
NConf = GetNConf(BaseName,NDigits);

Animation = NConf > 0;

if nargin == 0

   % Get the number of components:
   NC = NComp(filename);
   Opc = 2;
    while Opc == 2
        
        % Specify the resolution for the bodies:
        PR = input('Specify the particle resolution (default = 20): ');
        if isempty(PR)
            PR = 20;
        end
        
        % Specify the resolution for the walls:
        WR = input('Specify the wall resolution (default = 40): ');
        if isempty(WR)
            WR = 40;
        end
        
        % Specify the level of opacity of the walls:
        Opacity = input('Specify the opacity level of walls (default = 0.3): ');
        if isempty(Opacity)
            Opacity = 0.3;
        end
        
        if Animation
            
            % Specify the number of configurations to skip in each step:
            Skip = input('Specify the number of configurations to skip in each step (default = 0): ');
            if isempty(Skip)
                Skip = 0;
            end
            
            % Specify the number of frames per second:
            FPS = input('Specify the number of frames per second (default = 15): ' );
            if isempty(FPS)
                FPS = 15;
            end
            
        end
        
        % Set the default values:
        WallColor = [0.753 0.753 0.753];
        Color = rand(NC,3);
        
        % Open a dialog box for the user to choose the color of the walls:
        WallColor = uisetcolor(WallColor,'Choose a color for the walls:');
        if isequal(WallColor,0), return; end
        
        % Open dialog boxes for the user to choose component colors:
        for i = 1:NC
            C = uisetcolor(Color(i,:), ['Color for Component ' num2str(i)]);
            if isequal(C,0), return; end
            Color(i,:) = C;
        end
        
        % Specify if only half configuration must be ploted:
        Half = 2 - menu('Plot only half configuration?','Yes','No');
        
        % Specify if the objective is video or images:
        Video = 2 - menu('Which kind of result?','Video','Images');

        Fig = figure;
        Width = 800; Height = 600;
        Screen = get(0,'ScreenSize');
        ScreenWidth = Screen(3); ScreenHeight = Screen(4);
        Left = (ScreenWidth - Width)/2; Bottom = (ScreenHeight - Height)/2;
        set(Fig,'position',[Left Bottom Width Height]);
        set(Fig,'renderer','OpenGL');
        PlotConfig(filename,WallColor,Color,PR,WR,Opacity,3,Half);
        
        if Animation
            Opc = menu('You can now resize the window and/or rotate the axes.', 'Ok', 'Redefine','Cancel');
            
            if Opc == 0
                return
            elseif Opc == 1
                View = get(get(Fig,'CurrentAxes'),'View');
            elseif Opc == 3
                Opc = menu('Save figure?','Yes','No');
                if Opc == 1
                    SaveTiff( BaseName, Half );
                    title(['Saved as ' BaseName '.tif']);
                end
                return;
            else
                close(Fig);
            end
        else
            Opc = menu('Save figure?','Yes','No');
            if Opc == 1
                SaveTiff( BaseName, Half );
                title(['Saved as ' BaseName '.tif']);
            end
            return;
        end
    end
    
    if Video
        % Open a dialog box for the user to entre the name of the output file:
        [AVIfilename path] = uiputfile({'*.avi' 'Audio-Video Interleaved files (*.avi)'},'Save animation as...');
        if isequal(AVIfilename,0) | isequal(path,0)
            close(Fig);
            return;
        else
            AVIfilename = [path AVIfilename];
        end
    end

    PBC = [0 0 0];
    MaxNW = -1;

else
    if nargin < 15
        MaxNW = -1;
        if nargin < 14
            Video = 1;
            if nargin < 13
              disp('The number of arguments of the function Configuration must be 0 or at least 13.');
              return;
            end
        end
    end
end

if Video
    if Half
        AVIfilename = [BaseName '_Half.avi'];
    else
        AVIfilename = [BaseName '.avi'];
    end
    % Initialize AVI file:
    mov = avifile(AVIfilename,'compression','Indeo5','quality',100,'fps',FPS);
end

if nargin ~= 0
    Fig = figure;
    Screen = get(0,'ScreenSize');
    ScreenWidth = Screen(3); ScreenHeight = Screen(4);
    Left = (ScreenWidth - Width)/2; Bottom = (ScreenHeight - Height)/2;
    set(Fig,'position',[Left Bottom Width Height]);
    set(Fig,'renderer','OpenGL');
end

% Plot every configuration and get the resultant frame for the AVI file:
if Animation
    Index = [0:Skip+1:NConf];
    NIndex = length(Index);
    Aux = ['/' int2str(NIndex)];
    for i = 1 : NIndex
        
        Name = FileSeries(BaseName,Index(i),NDigits);
        
        clf;
        PlotConfig([Name '.cfg'],WallColor,Color,PR,WR,Opacity,View,Half,PBC,MaxNW);
        
        title(['Frame: ' int2str(i) Aux]);
        
        if Video 
            mov = addframe(mov,getframe);
        else
            SaveTiff( BaseName, Half );
        end
        
    end
    
    if Video, mov = close(mov); end
    
    close(Fig);
    
else
    PlotConfig(filename,WallColor,Color,PR,WR,Opacity,View,Half,PBC,MaxNW);
    if ~Video
        SaveTiff( BaseName, Half );
    end
end

%--------------------------------------------------------------------------
function [] = SaveTiff( BaseName, Half )

if Half
    print([BaseName '_Half.tif'],'-dtiff','-r300');
else
    print([BaseName '.tif'],'-dtiff','-r300');
end

%--------------------------------------------------------------------------
function BaseName = GetBaseName( Name )

Underscore = findstr('_',Name);
if isempty(Underscore)
  Ext = findstr('.',Name);
  if isempty(Ext)
      BaseName = Name;
  else
      BaseName = Name(1:Ext-1);
  end
else
    BaseName = Name(1:Underscore-1);
end

%--------------------------------------------------------------------------
function NDigits = GetNDigits(BaseName)

Continua = 1;
NDigits = 0;
while Continua
    NDigits = NDigits + 1;
    Continua = exist([FileSeries(BaseName, 0, NDigits ) '.cfg'],'file') == 0 & NDigits < 10;
end

%--------------------------------------------------------------------------
function FileName = FileSeries( Name, Index, NumberOfDigits )

for i = NumberOfDigits-1:-1:0
    D = Index;
    for j = NumberOfDigits-1:-1:i+1
        D = D - 10^j*Digit(j+1);
    end
    D = fix(D / 10^i);
    Digit(i+1) = D;
end

FileName = [Name '_'];
for i = NumberOfDigits:-1:1 
    FileName = [FileName int2str(Digit(i))];
end

%--------------------------------------------------------------------------
function NConf = GetNConf(BaseName,NDigits)

if exist([FileSeries(BaseName, 0, NDigits) '.cfg'],'file') ~= 0
    NConf = 0;
    while exist([FileSeries(BaseName, NConf + 1, NDigits) '.cfg'],'file') ~= 0
        NConf = NConf + 1;
    end
else
    NConf = -1;
end
