function  kmeans2(T, K, maxIter)
	# T is the data to be clustered, where rows represent datapoints and
	# columns the dimension of the data. C is the matrix with the centroids of
	# each cluster, here row represents the cluster and the columns are the
	# dimensions of the cluster centroid. W is the corresponding weights of
	# each cluster. maxIter is maximum number of iterations allowed. K is number of clusters.

	eps = 0.000005;
	(M, k) = size(T);
	
	#Kludgy sample rows without replacement
	rand_indicies = rand(1:size(T,1),100*k);
	rand_indicies = unique(rand_indicies); #order preserving duplicate deletion
	rand_indicies = rand_indicies[1:K];
	C = T[rand_indicies,:];

	dist_old = 1;
	iter = 0;
	delta = 1;
	w = zeros(1, K);
	

	while delta > eps && maxIter > iter
	    # D is distance and Q is its index in C, returns only the smallest
	    # euclidean distances.
	    p_dists = pdist(C, T);
	    p_dists=convert(Array{Float64,2},p_dists); #this is needed since otherwise findmin() fails on it.
	    (D,Q) = findmin(p_dists,1);
    	Q = map(x->ind2sub(size(p_dists), x)[1], Q); #Q is an ordered pair. I just want the location in the row.
    
	    dist = 1/(M*k)*sum(D.^2);
    
	    w = zeros(1, K);
	    C = zeros(K, k);
    
    	# Update the weights and recompute centroid of each cluster    
    	for i = 1:length(Q)
        	w[Q[i]] = w[Q[i]] + 1; #This is basically tallying up
	        C[Q[i],:] = C[Q[i],:] + T[i,:];
    	end

	    #Might be able to make this faster with something like (C'./w)'
    	for i = 1:length(w)
        	C[i,:] = C[i,:]/w[i];
	    end

    
	    w = w/sum(w);
    
    	iter = iter + 1;
    
	    delta = abs(dist - dist_old)/dist_old;
    	dist_old = dist;

	end
	# exception of zero probability cluster, turn the NANs to the zero vector
	C[find(isnan(C))] = 0;
	return (C,w)
end