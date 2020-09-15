# RDTool

### command line options

    -j, job type (MAPPING/CONSENSUS)
    
    -rdf_id, reaction id field (for RDF format)

### MAPPING arguments
	
	-i, input file
    
    -o, output directory
    
    -fmt, format (RDF/SMI), optional, default = RDF
    
    -t, optional; if specified, defines maximum number of seconds a substructures enumeration will run for
    
    -min, skip MIN algorithm, default = false
    
    -max, skip MAX algorithm, default = false
    
    -mixture, skip MIXTURE algorithm, default = false

*Example: java -jar release/aam-utils.jar -j MAPPING -i first_10.rdf -o first_10_results -rdf_id Reaction_ID*

### CONSENSUS arguments

    -i, list of input files in RDF format, separated with semicolon
    
    -o, output file in RDF format
    
    -ignore_tfc, ignore "total_fragment_changes" metric

*Example: java -jar release/aam-utils.jar -j CONSENSUS -i "first_10_results/MIN_reactions.rdf;first_10_results/MAX_reactions.rdf" -o first_10_results/consensus.rdf*

# dependencies
	- cdk: 2.3
	- commons-cli: 1.4
	- jgrapht-core: 1.1
