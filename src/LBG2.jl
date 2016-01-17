function LBG2(T, C, w)
	# T is the data to be clustered, where rows represent datapoints and
	# columns the dimension of the data. C is the matrix with the centroids of
	# each cluster, here row represents the cluster and the columns are the
	# dimensions of the cluster centroid. W is the corresponding weights of
	# each cluster.

	eps = 0.000005;
	(M, k) = size(T);
	N = length(w);

	# Squared-error distortion measure.
	#LGB2 uses pdist2 (which I renamed to pdist3) which I found at http://www.mathworks.com/matlabcentral/fileexchange/29004-feature-points-in-image--keypoint-extraction/content/FPS_in_image/FPS#20in#20image/Help#20Functions/SearchingMatches/pdist2.m
	dist = 1/(M*k)*sum(minimum(pdist(C, T),1));
	dist_old = 0;

	# Find index of the largest cluster.
	I = indmax(w);

	# Split it
	C = [C; zeros(1, k)]; #Add a row of zeros on the bottom
	C[I, :] = (1+eps)*C[I, :];
	C[N+1, :] = (1-eps)*C[I, :];
	N = N + 1;

	iter = 0;
	delta = 1;

	while delta > eps
    	# D is distance and Q is its index in C, returns only the smallest
	    # euclidean distances.
    	p_dists = pdist(C, T);
    	p_dists=convert(Array{Float64,2},p_dists); #this is needed since otherwise findmin() fails on it.
	    (D,Q) = findmin(p_dists,1);
    	Q = map(x->ind2sub(size(p_dists), x)[1], Q); #Q is an ordered pair. I just want the location in the row.
    
	    w = zeros(1, N);
    	C = zeros(N, k);
    
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
    	dist_old = dist;
	    dist = 1/(M*k)*sum(D.^2);
    	delta = (dist - dist_old)/dist_old;
	end
	return (C,w)
end