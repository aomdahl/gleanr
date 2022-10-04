filelist = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/trait_selections/infertility/fall_2022/sept_2022_list.txt" #get this from the YAML
shell.prefix("ml anaconda;")
#prep the summary stats for LDSC analysis.

#Split the calls and run separately for each trait so we get the sum stats the way we want it...
def getTraitNames(fl):
    traits = list()
    with open(fl, 'r') as istream:
        for l in istream:
            l = l.strip().split()
            traits.append(l[1])
    return(traits)

def getGWASNames(fl):
    gwas = list()
    import os
    with open(fl, 'r') as istream:
        for l in istream:
            l = l.strip().split()
            basename = os.path.basename(l[0])
            gwas.append(os.path.splitext(basename)[0])
    return(gwas)
trait_list=getTraitNames(filelist)

rule all:
	input:
		expand("ldsr_results/infertility/fall_2022/{trait_id}.log", trait_id = trait_list) 

rule mungeCommands:
    input: 
        study_list="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/trait_selections/infertility/fall_2022/sept_2022_list.txt", snp_ref="/data/abattle4/aomdahl1/reference_data/hm3_snps_ldsc_ukbb.tsv"

    output: "gwas_extracts/{identifier}/munge_sumstats.all.sh"
    params:
        "gwas_extracts/{identifier}/"
    shell:
        """
        conda activate std
        python src/munge_precursor.py --study_list {input.study_list} --merge_alleles {input.snp_ref} --output {params} --mod
        #we want the mod version...
        """


rule splitMungeCommands:
    input: "gwas_extracts/{identifier}/munge_sumstats.all.sh"
    output:
        expand("gwas_extracts/{{identifier}}/munge_calls/ldsc_line_runs.{runid}.sh", runid=trait_list)
    params:
        "gwas_extracts/{identifier}/munge_calls/ldsc_line_runs."
    run:
        import re
        ex_list = list()
        print(str(input[0]))
        with open(str(input[0]), 'r') as istream:
            for l in istream:
                txt = l.strip()
                try:
                    tag = re.search("--out [-\w\/\.]+\.([\w_]+) --", txt)
                    with open(params[0] + str(tag.group(1)) + ".sh", 'w') as ostream:
                        ostream.write(txt)
                        if "INTERMEDIATE." in txt:
                            intermediate_file=re.search("--sumstats (.+) --out", txt)
                            ostream.write("\n#rm " + intermediate_file.group(1))
                    ex_list.append(tag)
                except AttributeError:  
                    print("Error with REGEX- debug")



#use the expand to get through the entire list?
#th4 names basically come from the full file name less the extension. Not sure this is the way we want to do it...
#NEed to test this out...
rule ldscMunge:
    input: "gwas_extracts/{identifier}/munge_calls/ldsc_line_runs.{trait_id}.sh"
    output: 
        "gwas_extracts/{identifier}/{gwas_id}.{trait_id}.sumstats.gz"
    shell:
	    """
        conda activate py27
	    bash {input}
	    """

gwas_list = getGWASNames(filelist)
full_study_ids = [gwas_list[i] + "." + trait_list[i] for i in range(0,len(gwas_list))]



rule generateLDSCCommands:
    input: expand("gwas_extracts/{{identifier}}/{study_id}.sumstats.gz",study_id=full_study_ids)
    output:
        traitfile = "ldsr_results/{identifier}/pairwise.traits.txt",
    
    run:
        #Make a file with all the names for running
        import os
        wd=os.getcwd()
        with open(output.traitfile, 'w') as ostream:
            for i in input:
                ostream.write(wd + "/" + i + '\n') 
        #Generate all the commands. Puts each one in its own file...

rule pairwiseCommands:
    input:
        "ldsr_results/{identifier}/pairwise.traits.txt"
    output:
        expand("ldsr_results/{{identifier}}/{trait_id}_ldsc.run.sh", trait_id = trait_list)
    shell:
	    """
	     fp=`pwd`
	     bash src/all_pairwise_r2g.sh {input}  $fp/ldsr_results/{wildcards.identifier}/
	    """

#This actually runs the LDSC that is needed. Allows us to just run single ones as needed.
rule pairwiseLDSC:
    input:
        "ldsr_results/{identifier}/{trait_id}_ldsc.run.sh"
    output:
        "ldsr_results/{identifier}/{trait_id}.log"
    shell:
        """
		d=`pwd`
       		bash {input}
		cd $d
        """

#TODO: upgrade this so it can do liability scale
rule liabilityScaleHeritability: #not actually yet, need to do this
    input: 
    	"gwas_extracts/{identifier}/{study_id}.sumstats.gz"
    output:
        "ldsr_results/{identifier}/h2_ldsr/{study_id}.log"
    params:
        "ldsr_results/{identifier}/h2_ldsr/{study_id}"
    shell:
        """
        ml anaconda
        conda activate py27
        PTH=`pwd`
	cd /data/abattle4/aomdahl1/reference_data/ldsc_ref/
        python2  ~/.bin/ldsc/ldsc.py \
            --h2 $PTH/{input} \
            --ref-ld-chr eur_w_ld_chr/ \
            --w-ld-chr eur_w_ld_chr/ \
            --out $PTH/{params}
        #we'd like to be able to get things on the liability scale, huh....
        #--samp-prev 0.5,0.5 \
        #--pop-prev 0.01,0.01
        cd -
    """
