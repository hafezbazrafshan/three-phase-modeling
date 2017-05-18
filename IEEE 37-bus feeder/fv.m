function  [fv]=fv(vR,sLoad,Y2,vS)

vr=vR(1:3); 
vi=vR(4:6);

f=-diag(vr-j*vi)*Y2*(vr+j*vi)+diag(vr-j*vi)*Y2*vS-conj(sLoad);

fv=[real(f); imag(f)];
