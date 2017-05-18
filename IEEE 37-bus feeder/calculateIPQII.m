function fv= calculateIPQII(v, g, c, e, sL,iL,yL )


g(:,3)=0;  % since yL is already accounted for in the Ybus using getYLoadImpedance function
nvars=length(v); 
fv=zeros(nvars,1);
for i=1:nvars
    eMat=reshape(e(i,:,:),nvars,2);
    eVec1=eMat(:,1);
    eVec2=eMat(:,2);
    iL_PQ= c(i,:)*[ fPQ( eVec1.' * v, sL(i,1)); fPQ(eVec2.'*v, sL(i,2))];
    iL_I=c(i,:) *[ fI(eVec1.'*v, iL(i,1)); fI(eVec2.'*v, iL(i,2))]; 
    iL_Y=c(i,:)*[fY(eVec1.'*v, yL(i,1)); fY(eVec2.'*v, yL(i,2))];
    fv(i,1)=g(i,:)*[iL_PQ;iL_I;iL_Y];
end


end

