function [Bus,pu1,angle1,pu2,angle2,pu3,angle3]=importKerstingSolutions(filename)

if nargin < 1
    filename='powerflow.txt';
end
fileID=fopen(filename); 
tline = fgetl(fileID);
idx=0;  % counter for the power flows (flows on links)
ii=0;  % counter for nodes
while ischar(tline)
  %% ****Obtaining the node  
    match1=strfind(tline,'NODE:');
  
  %% **** Obtaining the voltages
  if ~isempty(match1)
      ii=ii+1;
  match2=strfind(tline,'VOLTS:'); 
  nodelabelString=tline(match1+5:match2-1);
  nodelabel=nodelabelString; 
  
      
      if length(strfind(nodelabelString, 'RG1'))>=1
          nodelabel='150r';
          
      elseif length(strfind(nodelabelString, 'RG2'))>=1
          nodelabel='9r';
          
      elseif length(strfind(nodelabelString, 'RG3'))>=1
          nodelabel='25r';
          
          
      elseif length(strfind(nodelabelString,'RG4'))>=1
          nodelabel='160r';
          
          
      elseif length(strfind(nodelabelString,'XF1'))>=1
          nodelabel='61s';
          
      end
      

  
  Bus{ii,1}=strtrim(nodelabel);
  pu1{ii,1}  =nanForEmptyCell(str2num(tline(match2+6:31)));
  angle1{ii,1}=nanForEmptyCell(str2num(tline(32:40)));
    pu2{ii,1}=nanForEmptyCell(str2num(tline(42:47)));
 angle2{ii,1}=nanForEmptyCell(str2num(tline(48:55)));
  pu3{ii,1}=nanForEmptyCell(str2num(tline(56:63)));
  angle3{ii,1}=nanForEmptyCell(str2num(tline(64:71)));
  end
  
  %% **** Obtaining from node power flows
  tline=fgetl(fileID); 
 
      
end
    

fclose(fileID);   