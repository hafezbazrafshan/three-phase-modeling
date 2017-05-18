
clear all;
clc;


fid=fopen('LineCodes2.dat','r'); 

sequences={'r1=','r0=','x1=','x0=','c1=','c0='};
sequencesStrings={'r1','r0','x1','x0','c1','c0'};
lineCnt=0;
lineDataSet={};

a=exp(sqrt(-1)*2*pi/3);
T=[1 1 1; 1 a.^2 a; 1 a a.^2];


tline = '';
while ~isequal(tline,-1)
    R=[];
    X=[];
    C=[];
sequenceImpedance=0;
    
    
    tline=fgetl(fid);
if strfind(tline,'[Linecode')>0
    
  
    
 tline = fgetl(fid);     

if strfind(tline,'name=')>0
    name=lower(strtrim(tline(strfind(tline,'name')+5:end)));
    tline=fgetl(fid);
end


if strfind(tline,'nphases=') >0
nphases=str2num(tline(strfind(tline,'nphases=')+8:end));
tline=fgetl(fid);
end

if strfind(tline,'units=') >0
units=strtrim(tline(strfind(tline,'units=')+6:end));
tline=fgetl(fid);
end


% if there are any sequences

if strfind(tline,'r1=') >0
r1=str2num(tline(strfind(tline,'r1=') +3:end));
r2=r1;
sequenceImpedance=1;
tline=fgetl(fid);
end

if strfind(tline,'r0=') >0
r0=str2num(tline(strfind(tline,'r0=') +3:end));
tline=fgetl(fid);
end

if strfind(tline,'x1=') >0
x1=str2num(tline(strfind(tline,'x1=') +3:end));
x2=x1;
tline=fgetl(fid);
end


if strfind(tline,'x0=') >0
x0=str2num(tline(strfind(tline,'x0=') +3:end));
tline=fgetl(fid);
end



if strfind(tline,'c1=') >0
c1=str2num(tline(strfind(tline,'c1=') +3:end));
c2=c1;
tline=fgetl(fid);
end


if strfind(tline,'c0=') >0
c0=str2num(tline(strfind(tline,'c0=') +3:end));
tline=fgetl(fid);
end

%************* RMATRIX
if isempty(strfind(name,'triplex'))
    R=zeros(nphases,nphases); 
    X=zeros(nphases,nphases); 
    C=zeros(nphases,nphases); 
    
if strfind(tline,'[Rmatrix]')
    tline=fgetl(fid); 
    R(1,1)=str2num(tline);
    tline=fgetl(fid);
    if isempty(strfind(tline,'[Xmatrix]'))
    cc=strsplit(tline); 
    R(1,2)=str2num(cc{1});
    R(2,1)=R(1,2);
    R(2,2)=str2num(cc{2});
    tline=fgetl(fid);
    end
    
    
    if isempty(strfind(tline,'[Xmatrix]'))
    cc=strsplit(tline); 
    R(1,3)=str2num(cc{1});
    R(3,1)=R(1,3);
    R(2,3)=str2num(cc{2});
    R(3,2)=R(2,3);
    R(3,3)=str2num(cc{3});
    tline=fgetl(fid);
    end
    
    
end


%************* XMATRIX
if strfind(tline,'[Xmatrix]')
    tline=fgetl(fid); 
    X(1,1)=str2num(tline);
    tline=fgetl(fid);
    if isempty(strfind(tline,'[Cmatrix]'))
    cc=strsplit(tline); 
    X(1,2)=str2num(cc{1});
    X(2,1)=X(1,2);
    X(2,2)=str2num(cc{2});
    tline=fgetl(fid);
    end
    
    
    if isempty(strfind(tline,'[Cmatrix]'))
    cc=strsplit(tline); 
    X(1,3)=str2num(cc{1});
    X(3,1)=X(1,3);
    X(2,3)=str2num(cc{2});
    X(3,2)=X(2,3);
    X(3,3)=str2num(cc{3});
    tline=fgetl(fid);
    end
    
end
    
    
    %************* CMATRIX
if strfind(tline,'[Cmatrix]')
    tline=fgetl(fid); 
    C(1,1)=str2num(tline);
    tline=fgetl(fid);
    if ~isempty(tline)
    cc=strsplit(tline); 
    C(1,2)=str2num(cc{1});
    C(2,1)=C(1,2);
    C(2,2)=str2num(cc{2});
    tline=fgetl(fid);
    end
    
    
    if ~isempty(tline)
    cc=strsplit(tline); 
    C(1,3)=str2num(cc{1});
    C(3,1)=C(1,3);
    C(2,3)=str2num(cc{2});
    C(3,2)=C(2,3);
    C(3,3)=str2num(cc{3});
    tline=fgetl(fid);
    end
    
    
    
end
else
end



if sequenceImpedance==0
Zseries=R+sqrt(-1)*X;
Yshunt=(sqrt(-1)*C*10^(-9));
else
Zseries= T* [r0+sqrt(-1)*x0 0 0; 0 r1+sqrt(-1)*x1 0; 0 0 r2+sqrt(-1)*x2]*inv(T);
Yshunt=sqrt(-1)*T * (1i*2*pi*60*[c0 0 0; 0 c1 0; 0 0 c2])*inv(T)*10^(-9);

if nphases==1
    Zseries=Zseries(1,1); 
    Yshunt=Yshunt(1,1); 
end
end

lineCnt=lineCnt+1;
lineDataSet{lineCnt,1}=name;
lineDataSet{lineCnt,2}=nphases;
lineDataSet{lineCnt,3}=Zseries;
lineDataSet{lineCnt,4}=Yshunt;
end
end


        
fclose(fid);

save('lineDataSet2','lineDataSet'); 

 
