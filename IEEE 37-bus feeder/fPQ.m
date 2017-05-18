
function out = fPQ( v,s)
% FPQ   Constant PQ Load i=fPQ(v,s) outputs the single phase current
% i, given the single phase inputs v, and the single phase apparent power
% consumption s

if abs(v) >0 
out=conj(s./v) ;

else 
    out=0;
end





end

