# ==============================================================================
# ARKQuikr.jl
#
# Authors: David Koslicki (david.koslicki@math.oregonstate.edu)
#
# Takes a fasta file as input, performs the ARK algorithm on it,
# and then performs the Quikr Algorithm on it.
# ==============================================================================

###################
#To do:
#1. Add the SEK database option
#####################

using ArgParse
using HDF5
include("lsqnonneg.jl")
include("LBG2.jl")
include("pdist.jl")
include("ConvertToCAMIOutput.jl")
include("kmeans2.jl")



#Parse arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input_file", "-i"
			help = "Fasta file to perform the reconstruction on"
        "--lambda", "-l"
			help = "Lambda parameter. The lower it is, the sparser the reconstruction. The higher it is, the more closely the kmer counts will be fit. Default is 10,000"
			default = 10000
		"--output_file", "-o"
			help = "Output text file"
		"--kmer_counts_per_sequence_path", "-k"
			help = "Full path to the kmer_counts_per_sequence file obtained from https://github.com/mutantturkey/dna-utils"
			default = "kmer_counts_per_sequence"
		"--training_database", "-t"
			help = "Training database to use. Choices are: SEK, Quikr. SEK is more accurate, but slower, whereas Quikr is faster, but less accurate. Default is Quikr"
			default = "Quikr"
		"--number_of_clusters", "-n"
			help = "The number of clusters to use (at max). Default = 10"
			default = 10
		"--clustering_type", "-t"
			help = "Clustering type: Deterministic (slow), Random (fast). Default is Random"
			default = "Random"
			
    end
    return parse_args(s)
end


#Parse the args
parsed_args = parse_commandline()
input_file = parsed_args["input_file"]
lambda = int(parsed_args["lambda"])
output_file = parsed_args["output_file"]
kmer_counts_per_sequence_path = parsed_args["kmer_counts_per_sequence_path"]
training_database = parsed_args["training_database"]
number_of_clusters = int(parsed_args["number_of_clusters"])
clustering_type = parsed_args["clustering_type"]

k=6;
#Form the 6mer counts
counts_per_sequence_mult = map(x->int(split(strip(x))),readlines(`$kmer_counts_per_sequence_path -i $input_file -k $k -c`));
counts_per_sequence = [counts_per_sequence_mult[i][j] for i=1:length(counts_per_sequence_mult),j=1:4^k];

#Normalize the rows
counts_per_sequence_sums = sum(counts_per_sequence,2);
for i=1:size(counts_per_sequence,1)
	counts_per_sequence[i,:] = counts_per_sequence[i,:]/counts_per_sequence_sums[i];
end
#Type conversion
counts_per_sequence=convert(Array{Float64,2},counts_per_sequence);

#Now run the algorithm
if training_file == "Quikr" #Using the quikr database
	#Read in the training database
	A = h5read("../../data/trainset7_112011N6C.h5","/data");

	#Form the Aaux
	Aaux = [ones(1,size(A,2)); lambda*A];

	if  clustering_type == "Random"
		#Now for ARK
		# Initialization
		NoOfClusters_Quikr = number_of_clusters;
		
		(C_ARK_Quikr, ClusterProbability) = kmeans2(counts_per_sequence, NoOfClusters_Quikr, 1000); #Max of 1000 iterates of the clustering algorithm
		Mu_ARK_Quikr = C_ARK_Quikr';  #Cluster mean vectors
	    result_ARK_Quikr = zeros(1,size(A,2));
    	
	    #Perform Quikr on each cluster
    	for i=1:NoOfClusters_Quikr
        	s = [0; lambda*Mu_ARK_Quikr[:,i]];
	        tmp_ARK_Quikr = lsqnonneg(Aaux, s)';
    	    result_ARK_Quikr = result_ARK_Quikr + ClusterProbability[i]*tmp_ARK_Quikr; # This is the linear additive composition estimation
	    end
    
    	#Normalize the solution
	    result_ARK_Quikr = result_ARK_Quikr/sum(result_ARK_Quikr);
	
	elseif clustering_type == "Deterministic"
		#Initialization
		eta=0.005;
		Composition_ARK_Quikr  = zeros(1,size(A,2));
		ChangeInComposition_ARK_Quikr = 1;
		NoOfClusters_Quikr = 0;
		MaxNoOfClusters = number_of_clusters;
		
		while (ChangeInComposition_ARK_Quikr  > eta) && (NoOfClusters_Quikr < MaxNoOfClusters)  # (stopping criteria for LBG based clustering)
		    
		    #Perform the clustering
		    if NoOfClusters_Quikr == 0
    		    C_ARK_Quikr = mean(counts_per_sequence,1);
    		    ClusterProbability = [1.0];
	    	else
    	    	(C_ARK_Quikr, ClusterProbability) = LBG2(counts_per_sequence, C_ARK_Quikr, ClusterProbability);  # LBG algorithm increases the number of clusters as output from the input no of clusters by one 
		    end    
    		NoOfClusters_Quikr = length(ClusterProbability);    
    		
	    	# After clustering, Quikr is used for each cluster
	    	Mu_ARK_Quikr = C_ARK_Quikr';  # Cluster mean vectors)
		    result_ARK_Quikr = zeros(1,size(Aaux,2));    	
		    
		    #Perform Quikr on each cluster
    		for i=1:NoOfClusters_Quikr
	        	s = [0; lambda*Mu_ARK_Quikr[:,i]]; 
	        	tmp_ARK_Quikr = lsqnonneg(Aaux, s)';
    		    result_ARK_Quikr = result_ARK_Quikr + ClusterProbability[i]*tmp_ARK_Quikr; # This is the linear additive composition estimation
	    	end
	    	
	    	#Record the change in the composition
		    if NoOfClusters_Quikr > 1
        		ChangeInComposition_ARK_Quikr  = norm(Composition_ARK_Quikr[end,:] - result_ARK_Quikr, 1);
	    	end
	    	
	    	#Save the composition results
			Composition_ARK_Quikr  = [Composition_ARK_Quikr; result_ARK_Quikr]; 
		end
		result_ARK_Quikr = Composition_ARK_Quikr[end,:];
	else
		error("Invalid clustering_type. Choose one of: Deterministic, Random")
	end
    
    	#write the solution to file
		output_level = 0; #Since we don't have hypothetical organisms
		ConvertToCAMIOutput(result_ARK_Quikr, "../../data/trainset7_taxonomy.txt", output_level, output_file)
	
elseif training_file == "SEK" #using the split Quikr database (known as the SEK database)

	#Read in the training database
	A = h5read("../../data/trainset7_112011_allseqslongerthan700-SEKTrainingMatrix-bitShift100-windowLength400-N6C.h5","/data");
	
	#Read in the block matrix to return to the trainset7_SEK basis
	blockMatrix = h5read("../../data/trainset7_112011_allseqslongerthan700-SEKTrainingMatrix-bitShift100-windowLength400-blockMatrix.h5","/data");

	#Form the Aaux
	Aaux = [ones(1,size(A,2)); lambda*A];
	yaux = [0;lambda*counts];

	#Perform the reconstruction
	x = lsqnonneg(Aaux, yaux);
	
	#Return to the trainset7_SEK basis
	x = blockMatrix*x;

	#Normalize the output
	x = x/sum(x);

	#Write the output to file
	output_level = 0; #Since we don't have hypothetical organisms
	ConvertToCAMIOutput(x, "../../data/trainset7_SEK_taxonomy.txt", output_level, output_file)
end










#How I saved the training data
#fid = h5open("trainset7_112011N6C.h5","w");
#fid["/data", "chunk", (512,10046), "compress", 7] = Atemp;
#close(fid)

