function lsqnonneg(C::Matrix, d::Vector, tol::Real=-1, itmax_factor::Real=3)
	#set the tolerance
	(m,n) = size(C);
	tol = (tol == -1) ? tol = 10*eps()*norm(C,1)*(maximum(size(C))+1) : tol
	itmax = itmax_factor*n;
	
	# Initialize vector of n zeros and Infs (to be used later)
	wz = zeros(n,1);
	
	# Initialize set of non-active columns to null
	P = falses(n,1);
	
	# Initialize set of active columns to all and the initial point to zeros
	Z = trues(n,1);
	x = zeros(n,1);
    Ctrans=C';
    dtemp=d[:];
    BLAS.gemv!('N',-1.,C,x[:],1.,dtemp); #resid = d - C*x;
    w=BLAS.gemv('N',1.,Ctrans, dtemp); #w = Ctrans*resid;
	
	# Set up iteration criterion
	outeriter = 0;
	iter = 0;
	exitflag = 1;
	
	# Outer loop to put variables into set to hold positive coefficients
	while any(Z) & any(w[Z[:]] .> tol)
		#print("On iteration $(outeriter)\n")

		outeriter = outeriter + 1;
		
		# Reset intermediate solution z
		z = zeros(n,1);
		
		# Create wz, a Lagrange multiplier vector of variables in the zero set.
		# wz must have the same size as w to preserve the correct indices, so
		# set multipliers to -Inf for variables outside of the zero set.
		wz[P] = -Inf;
		wz[Z] = w[Z[:]];

		# Find variable with largest Lagrange multiplier
		t = indmax(wz);
		
		# Move variable t from zero set to positive set
		P[t] = true;
		Z[t] = false;
		
		# Compute intermediate solution using only variables in positive set
		z[P] = C[:,find(P)]\d;
		
		#inner loop to remove elements from the positive set which no longer belong
		while any(z[P] .<= 0)
			#print("entering inner loop\n")
			iter = iter + 1;
			if iter > itmax
				print("lsqnonneg:IterationCountExceeded");
		              exitflag = 0;
				iterations = outeriter;
		              resnorm = sum(resid.*resid);
				x = z;
				lambda = w;
				return x
			end
			
			# Find indices where intermediate solution z is approximately negative
			Q = [(z .<= 0) & P];

			# Choose new x subject to keeping new x nonnegative
			alpha = minimum(x[Q]./(x[Q] - z[Q]));
			x = x + alpha*(z - x);

			# Reset Z and P given intermediate values of x
			Z = [((abs(x) .< tol) & P) | Z ];
			P = ~Z;
			z = zeros(n,1);           # Reset z
			z[P] = C[:,find(P)]\d;      # Re-solve for z
		end

		x = z;

        dtemp=d[:];
        BLAS.gemv!('N',-1.,C,x[:],1.,dtemp); #resid = d - C*x;
        w=BLAS.gemv('N',1.,Ctrans, dtemp); #w = Ctrans*resid;

	end
	return x
end

