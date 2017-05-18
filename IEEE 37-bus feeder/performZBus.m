function [ network ] =performZBus( network, maxIt)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

if nargin <2
    maxIt=10; 
end

v2struct(network); 
v2struct(noLoadQuantities);
v2struct(loadQuantities);
availableBusIndices=network.availableBusIndices;

vPr=repmat(v0,N,1);   
 vPr=vPr(availableBusIndices(1:end-3));

success=0;
maxIt=10;
 

err=inf+zeros(maxIt+1,1);


vIterations=zeros(size(vPr,1),maxIt+1);
vIterations(:,1)=vPr;

fv= calculateIPQII(vPr, gMat, cMat, ePage, sL_load,iL_load,yL_load );

 err(1)=sum(abs((Y+YL)*vPr+fv+Y_NS*v0));
str=['Iteration No. ', num2str(0),' Error is ', num2str(err(1)), '\n'];
fprintf(str);
 for it=1:maxIt
     
      if err(it) <1e-5
     itSuccess=it-1;
     fprintf('Convergence \n');
     success=1;
     break;
      end 
      
      
      


vNew= inv(Y+YL)* (-fv-Y_NS*v0);
vIterations(:,it+1)=vNew;
 vPr=vNew;
 fv= calculateIPQII(vPr, gMat, cMat, ePage, sL_load,iL_load,yL_load );
 err(it+1)=sum(abs((Y+YL)*vNew+fv+Y_NS*v0));
 str=['Iteration No. ', num2str(it),' Error is ', num2str(err(it+1)), '\n'];
 fprintf(str);
 end
 
 
 vsol=vNew;
 
vsol=[vsol;v0];
v=zeros(3*(N+1),1);

v(availableBusIndices,1)=vsol;
v(missingBusIndices,1)=NaN;




network.ZBusResults=v2struct(v, vsol, success, err); 




end

