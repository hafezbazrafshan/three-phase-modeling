function out =fI( v,i )
% FI   Constant Current Load i=fI(v,i) outputs the single phase current
% i, given the single phase inputs v, and the nominal current i

% if abs(v) >0
% out= i*v./abs(v);  % constant power factor model
% else 
% out=0;
% end
out=i;   % constant current phasor model

end

