function [N] = NComp(filename)

file = fopen(filename,'r');

ncc = 2;
huge = 2^(8*ncc-1)-2;
Complete = 72;

ib = 4;
rb = 8;

ncc = ['int' num2str(8*ncc)];
ib = ['int' num2str(8*ib)];
rb = ['float' num2str(8*rb)];

Signature = fread(file,1,'int32');
if ~isequal(Signature,1129137476), disp('Error: Invalid configuration file.'); return; end
Version = fread(file,1,'int8');
dim = fread(file,1,'int8');
Form = fread(file,1, 'int8');

Full = Form == Complete;

LengthScale = fread(file,1,rb);
TimeScale   = fread(file,1,rb);
if Full, MassScale = fread(file,1,rb); end
Time = fread(file,1,rb);
Lbox = fread(file,dim,rb);

N = fread(file,1,ib);

file = fclose(file);
