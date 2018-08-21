function [cost,assignments] = glomeruliAssignmentCostFunction(params)
    assignments=params.assignments; % vector of cluster assignments
    odorDist=params.odorDist;
    physDist=params.physDist;
    shapePriorNorm=params.shapePriorNorm;
    
    % cost function associated with assignments
    compositeDist = log10(odorDist) + log10(physDist) - log10(shapePriorNorm);
    
    cost=0;
    for j=1:size(compositeDist,2)
       cost=cost+compositeDist(j,assignments(j)); 
    end
end