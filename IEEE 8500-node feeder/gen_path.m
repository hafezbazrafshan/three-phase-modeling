function path = gen_path(n, inc_mat, sNodes, rNodes)

tempPath=[n];

children=rNodes(inc_mat(n,:)==1);

for idx=1:length(children)
    childIdx=children(idx);
    tempPath=[tempPath, gen_path(childIdx,inc_mat, sNodes, rNodes)];
end
 
path=tempPath;

        
