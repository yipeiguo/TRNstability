% This is a function for creating a random directed acyclic graph given the
% number of nodes and density of the interaction.
% Here we assume the the second last and last nodes represent RNAP and
% ribosome respectively and they don't regulate gene copy number.
% Here, we construct a DAG which by definition cannot have self-loops. But
% we know that transcription factors may affect themselves so these will be
% added in the main script.
% Note that the maximum number of connections is N*(N-1)/2.

function g = CreateRandomDAG(N,density,maxiter)

numInt = ceil(N*(N-3)*density);
g = sparse([],[],true,N,N);
posedges = 1:N*N;
posedges(1:N+1:end) = [];
posedges(mod(posedges,N)==N-1) = [];
posedges(mod(posedges,N)==0) = []; 
remainingedges = posedges;
iter = 0;
numedges = 0;
while nnz(g) < numInt
  edge = randsample(remainingedges,1); % get a random edge
  g(edge) = true;
  g(edge) = graphisdag(g);
  if g(edge) == true
      numedges = numedges + 1;
      remainingedges = remainingedges(remainingedges~=edge);
      %fprintf('number of edges: %d \n', numedges);
  end
  iter = iter+1;
  if iter > maxiter
      warning('maximum number of iterations exceeded!');
      break
  end
end

end