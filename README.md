# iterative-merlin

You can run the program as:

> ./MERLINNCA -d data.txt -r all_tfs.txt -g all_tgts.txt -p prior_M.txt -l 0.005 -t 1 -n 10 -o out
> 
> -d	is the expression matrix
> -r  is the list of TFs (we only consider these)
> -g  is the list of targets (we only consider these)
> -p  is the prior network
> -l  is the lambda for NCA
> -c  is the number of CV folds (instead of fixed lambda)
> -t  is the number of threads to use. 
>     If we use more than 1, the programs will break the target gene set into that many sets and run it on each set separately.
> -n  is the number of NCA replicates.
> -o  is output directory.

Sample inputs are mentioned at the end.

The method will infer TFA from prior, then infer a network from the TFA, and repeat 5 times.
The output directory will look like this:

> out/
> out/merged_${k}.txt		where k is the iteration index (from 0 to 5).
> 						If we have multiple threads, these will be concatenation of individual networks.
> 
> out/$k/					This will be the outputs for iteration k
> 
> out/$k/adj.txt			These are the NCA outputs
> out/$k/tfa.txt
> 
> out/$k/$i/				These will be the MERLIN outputs, for individual threads. 
> 						If we have 1 thread, there will be only a "0" directory.
> 
> 
> out/$k/$i/fold0/prediction_k300.txt			These are the outputs of MERLIN.
> out/$k/$i/fold0/modules.txt  


We have some sample input files:

> sample_data/:
> 	tfs.txt						All TFs in MacIsaac
> 	tfs_20.txt					20 TFs
> 	tgts.txt					All targets in MacIsaac
> 	tgts_100.txt				100 of them
> 	tgts_150.txt				150 of them
> 	clusters_tgts.txt			All gene in MacIsaac with cluster assignment
> 	clusters_tgts_600.txt		Cluster assignment for 600 genes
> 	exp.txt						Expression matrix (NatVar)
> 	exp_${k}.txt				Expression for kth sub-sample
> 	mac2.bszc.txt				MacIsaak network
> 	prior_M.txt					Motif prior
> 
> new_in/:
> 	all_tfs.txt					All TFs in NAR paper
> 	all_tgts.txt				All targets

