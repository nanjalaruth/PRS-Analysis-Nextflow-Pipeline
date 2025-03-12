nextflow.enable.dsl=2

process download_score_files {
    tag "Downloading score file ${blood_trait}_${pgs_id}"
    
    input:
        tuple val(blood_trait), val(pgs_id)

    output:
        tuple val(blood_trait), val(pgs_id), path("${pgs_id}.txt")

    script:
        """
        curl -O https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/${pgs_id}/ScoringFiles/${pgs_id}.txt.gz
        gunzip ${pgs_id}.txt.gz
        """
        
}

process download_meta_files {
    tag "Downloading meta file ${blood_trait}_${pgs_id}"
    
    input:
        tuple val(blood_trait), val(pgs_id)

    output:
        tuple val(blood_trait), val(pgs_id), path("${pgs_id}_metadata_scores.csv")

    script:
        """
        curl -O https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/${pgs_id}/Metadata/${pgs_id}_metadata_scores.csv
        """
}

process modify_meta_file {
    tag "Modifying meta file for ${blood_trait}_${pgs_id}"
    
    input:
        tuple val(blood_trait), val(pgs_id), path(meta_files)

    output:
        tuple val(blood_trait), val(pgs_id), path("${blood_trait}_${pgs_id}_meta_modified.txt")

    script:
        template "modify_meta_file.R"
}

process modify_score_file {
    tag "Modifying score file for ${blood_trait}_${pgs_id}"
    
    input:
        tuple val(blood_trait), val(pgs_id), path(score_files)

    output:
        tuple val(blood_trait), val(pgs_id), path("${blood_trait}_${pgs_id}_edit.txt")
    script:
        """
        #remove header information
        grep -v '^#' ${score_files} > ${blood_trait}_${pgs_id}_edit.txt
        """
}

process modify_score_file_2 {
    tag "Modifying score file_2 for ${blood_trait}_${pgs_id}"
    
    input:
        tuple val(blood_trait), val(pgs_id), path(score_files), path(hg19ref), path(hg38ref), path(reference_info)

    output:
        tuple val(blood_trait), val(pgs_id), path("${blood_trait}_${pgs_id}_edit_2.txt")

    script:   
        """
        #!/bin/bash
        # Get the genome build for the given PGS ID
        genome_build=\$(awk -F'\\t' -v id="\${pgs_id}" '\$2 == id {print \$3}' "${reference_info}")

        # Select the appropriate reference file based on genome build
        if [[ "\${genome_build}" == "hg38" || "\${genome_build}" == "GRCh38" ]]; then
            ref="${hg38ref}"
        elif [[ "\${genome_build}" == "hg19" || "\${genome_build}" == "GRCh37" ]]; then
            ref="${hg19ref}"
        else
            echo "Warning: Genome build not found for PGS ID \${pgs_id}, defaulting to hg19."
            ref="${hg19ref}"
        fi

        # Debugging: Check if the reference file is correctly assigned
        echo "Using reference file: \${ref}"

        # Check if the columns chr_name and chr_position exist in the input file
        if awk -F'\\t' 'NR==1 {c1=0; c2=0; for(i=1;i<=NF;i++){if(\$i=="chr_name"){c1=1}if(\$i=="chr_position"){c2=1}}} END{exit !(c1 && c2)}' "${score_files}"; then
            # If the file already contains chr_name and chr_position, append the SNP column
            awk -F'\\t' 'NR==1 {
                for(i=1; i<=NF; i++) { 
                    if (\$i == "chr_name") { chr_name_index=i } 
                    else if (\$i == "chr_position") { chr_position_index=i } 
                }
                print \$0, "\\tSNP";
                next
            } 
            { print \$0 "\\t" \$(chr_name_index) ":" \$(chr_position_index) }' "${score_files}" > "${blood_trait}_${pgs_id}_edit_2.txt"
        else
            # Extract headers separately
            header=\$(head -n 1 "${score_files}")

            # Sort files (excluding headers)
            tail -n +2 "${score_files}" | sort -k1,1 > "${score_files}_sorted"
            tail -n +2 "\${ref}" | sort -k3,3 > "\${ref}_sorted"

            # Perform join operation
            join -1 1 -2 3 -t \$'\\t' "${score_files}_sorted" "\${ref}_sorted" | \\
            awk '{print \$0 "\\t" \$(NF-1) ":" \$NF}' > "${blood_trait}_${pgs_id}_edit_2_body.txt"

            # Add back the header row
            echo -e "\${header}\\tchr_name\\tchr_position\\tSNP" | cat - "${blood_trait}_${pgs_id}_edit_2_body.txt" > "${blood_trait}_${pgs_id}_edit_2.txt"

            # Clean up temporary files
            # rm "${score_files}_sorted" "\${ref}_sorted" "${blood_trait}_${pgs_id}_edit_2_body.txt"
        fi
        #rm ${hg38ref} ${hg19ref}
        """
}

process liftover {
    tag "Liftover ${blood_trait}_${pgs_id}"
    
    input:
        tuple val(blood_trait), val(pgs_id), path(score_files), path(hg19ref), path(hg38ref), path(reference_info), val(tgt_genome_build), path(chain_hg19_to_hg38), path(chain_hg38_to_hg19) 

    output:
        tuple val(blood_trait), val(pgs_id), path("${blood_trait}_${pgs_id}_edit_3.txt")

    script:
        """
        #!/bin/bash
        set -e  # Exit on error

        # Reference genome build for the given PGS ID
        reference_build=\$(awk -F'\\t' -v id=${pgs_id} '\$2 == id {print \$3}' ${reference_info})
        target_build=${tgt_genome_build}

        if [[ "\$target_build" == "hg19" || "\$target_build" == "GRCh37" ]]; then
            echo "Target build: hg19"
            
            if [[ "\$reference_build" == "hg19" || "\$reference_build" == "GRCh37" ]]; then
                echo "Reference data is already in hg19. No conversion needed."
                cp ${score_files} ${blood_trait}_${pgs_id}_edit_3.txt
                exit 0
            else
                echo "Converting from hg38 to hg19..."
                chain_file=${chain_hg38_to_hg19}
            fi

        elif [[ "\$target_build" == "hg38" || "\$target_build" == "GRCh38" ]]; then
            echo "Target build: hg38"
            
            if [[ "\$reference_build" == "hg38" || "\$reference_build" == "GRCh38" ]]; then
                echo "Reference data is already in hg38. No conversion needed."
                cp ${score_files} ${blood_trait}_${pgs_id}_edit_3.txt
                exit 0
            else
                echo "Converting from hg19 to hg38..."
                chain_file=${chain_hg19_to_hg38}
            fi
        else
            echo "Invalid target genome build specified!"
            exit 1
        fi

        # Prepare the SNPs for liftOver (convert to BED format)
        awk -F'\\t' '
        NR==1 {
            for (i=1; i<=NF; i++) {
                if (\$i == "chr_name") chr_col = i;
                if (\$i == "chr_position") pos_col = i;
                if (\$i == "SNP") snp_col = i;
            }
            next
        }
        NR>1 {
            print "chr"\$chr_col, \$pos_col-1, \$pos_col, \$snp_col
        }' OFS="\\t" ${score_files} > snps_to_convert.bed

        # Perform liftOver
        if [ -z "\$chain_file" ]; then
            echo "Error: Chain file not set. Exiting."
            exit 1
        fi

        liftOver snps_to_convert.bed \${chain_file} ${blood_trait}_${pgs_id}_converted_snps.bed unmapped_snps.txt
        if [[ \$? -ne 0 ]]; then
            echo "Error: liftOver conversion failed."
            exit 1
        fi

        # Add header to the converted SNPs file
        echo -e "lifted_chr_name\\tstart_position\\tend_position\\tSNP" | cat - ${blood_trait}_${pgs_id}_converted_snps.bed > ${blood_trait}_${pgs_id}_edit_3.txt
        
        """
}

process modify_score_file_3 {
    tag "Modifying score file for ${blood_trait}_${pgs_id}"
    
    input:
        tuple val(blood_trait), val(pgs_id), path(edit_2), path(edit_3)

    output:
        tuple val(blood_trait), val(pgs_id), path("${blood_trait}_${pgs_id}_modified.txt")

    script:
        template "modify_score_file.R"
}


process compute_pgs_scores {
    tag "Computing PGS scores in PLINK for ${blood_trait}_${pgs_id}"
    publishDir "${params.outdir}/${dataset}/${blood_trait}", mode: 'copy', overwrite: true
    
    input:
        //tuple val(blood_trait), val(pgs_id), path(modified_scoring_file), val(dataset), path(bgen), path(index), path(sample)
        tuple val(blood_trait), val(pgs_id), path(modified_scoring_file), val(dataset), path(bed), path(bim), path(fam)

    output:
        tuple val(dataset), val(blood_trait), val(pgs_id), path("${dataset}_${blood_trait}_${pgs_id}_prsval*")

    script:
        base = bed.baseName
        """
        plink --bfile ${base} \\
        --score ${modified_scoring_file} 1 2 3 sum \\
        --out ${dataset}_${blood_trait}_${pgs_id}_prsval
        """
}


process modify_pgs_scores {
    tag "Modifying score file from PGS catalogue"
    
    input:
        tuple val(dataset), val(blood_trait), val(pgs_id), path(pgs_score)

    output:
        tuple val(dataset), val(blood_trait), val(pgs_id), path("${dataset}_${blood_trait}_${pgs_id}_pgs_scores.txt")

    script:
        """
        #extract FID and SCORE columns
        awk '{print \$1 "\\t" \$6}' ${pgs_score} > ${dataset}_${blood_trait}_${pgs_id}_pgs.txt
        awk -v OFS='\\t' '{if (NR==1) {\$2="${pgs_id}_SCORE"}; print}' ${dataset}_${blood_trait}_${pgs_id}_pgs.txt > ${dataset}_${blood_trait}_${pgs_id}_pgs_scores.txt
        """
        
}

process conc_scores {

    publishDir "${params.outdir}/${blood_trait}", mode: 'copy', overwrite: true
    
    tag "Concatenating pgs scores"
    
    input:
        tuple val(dataset), val(blood_trait), path(score)

    output:
        tuple val(dataset), val(blood_trait), path("${dataset}_${blood_trait}_pgs_scores.txt")

    script:
        template "merge_scores.R"
}

process intergrate_scores {
    tag "Modifying score file from PGS catalogue"
    publishDir "${params.outdir}/${blood_trait}", mode: 'copy', overwrite: true
    
    input:
        tuple val(blood_trait), path(score_files), path(pheno), path(cov)

    output:
        tuple val(blood_trait), path("${blood_trait}_protsignature.txt"),
         path("${blood_trait}_correlation.pdf"), path("${blood_trait}_pred_score.txt"),
         path("${blood_trait}_ug_prs_pcs.txt"), path("${blood_trait}_predicted.txt")

    script:
        template "intergrate_scores.R"
}

process association_analysis {
    tag "Association analysis"
    publishDir "${params.outdir}/${blood_trait}", mode: 'copy', overwrite: true
    
    input:
        tuple val(blood_trait), path(ug_prs_pcs), path(lipid_trait)

    output:
        tuple val(blood_trait), path("${blood_trait}_forest_plot.pdf"), 
            path("${blood_trait}_results.csv")

    script:
        template "prediction.R"
}

workflow{
     //step 1
     // Download score files from PGS catalogue
    input = Channel.fromList(params.pheno)
    //input.view()
    download_score_files(input)

    download_meta_files(input)
    //modify metadata file
    modify_meta_file(download_meta_files.out)

    //step 2
    // Modify score file
    modify_ch = download_score_files.out
        //.map{ btrait, pgsid, score_file -> 
          //          return [btrait, pgsid, score_file, pheno] 
            //    }
    //modify_ch.view()
    modify_score_file(modify_ch)

    //step 2.2
    hg19ref = Channel.fromPath(params.ref19)
    hg38ref = Channel.fromPath(params.ref38)
    input = modify_score_file.out
    final_input = input
        .combine(hg19ref)
        .combine(hg38ref)
        .combine(modify_meta_file.out, by:[0,1])
    //input.view()
    modify_score_file_2(final_input)

    // liftover hg38 to hg19
    tgt_genome_build = Channel.value(params.target_genome_build)
    chain_hg19_to_hg38 = Channel.fromPath(params.chain_hg19_to_hg38)
    chain_hg38_to_hg19 = Channel.fromPath(params.chain_hg38_to_hg19)
    input = modify_score_file_2.out
        .combine(hg19ref)
        .combine(hg38ref)
        .combine(modify_meta_file.out, by:[0,1])
        .combine(tgt_genome_build)
        .combine(chain_hg19_to_hg38)
        .combine(chain_hg38_to_hg19)
    //input.view()
    liftover(input)

    //step 2.3
    modify_3_ch = modify_score_file_2.out
    .combine(liftover.out, by:[0,1])
    //modify_3_ch.view()
    modify_score_file_3(modify_3_ch)

    //step 3
    // calculate pgs score
    modify_score_out = modify_score_file_3.out
    plink_ch = Channel.fromList(params.plink_file)
    pgs_ch = modify_score_out.combine(plink_ch)
    //pgs_ch.view()
    compute_pgs_scores(pgs_ch)

    //// for batch computing
    ////modify_score_out = modify_score_file_3.out
    ////plink_ch = Channel
    ////.from(params.batches)
    ////.map { batch ->
        //// Construct file paths based on the naming format
        ////def bed_file = "${params.path}/${params.prefix}${batch}${params.suffix}.bed"
        ////def bim_file = "${params.path}/${params.prefix}${batch}${params.suffix}.bim"
        ////def fam_file = "${params.path}/${params.prefix}${batch}${params.suffix}.fam"
        
        //// Return the batch and associated file paths as a list
        ////[batch, bed_file, bim_file, fam_file]
    ////}
    ////pgs_ch = modify_score_out.combine(plink_ch)
    ////pgs_ch.view()
    ////compute_pgs_scores(pgs_ch)
   

    //step 3.2 modify_pgs_scores
    in = compute_pgs_scores.out
     .map{dset, btrait, pgsid, pgs_score -> [dset, btrait, pgsid, pgs_score[2]]}
     ////in.view()
    modify_pgs_scores(in)

    //step 3.3 conc scores
    input = modify_pgs_scores.out
      //  //.map{dset, btrait, pgsid, pgs_score -> [btrait, pgs_score]}
        //// .groupTuple()
        .map{dset, btrait, pgsid, pgs_score -> [dset, btrait, pgs_score]}
        .groupTuple(by: [0, 1])
       
    //input.view()
    conc_scores(input)

    ////step 4 Intergrate scores
    ////scores = conc_scores.out
    ////scores.view()
    ////pheno = Channel.fromList(params.phenotype)
    //cov = Channel.fromPath(params.covariate)
    //cov_pheno = pheno
      //  .combine(cov)    
    //test_in = scores
      //  .combine(cov_pheno, by:0)
    //intergrate_scores(test_in)

    //step 5 Predictions
    //trait = Channel.fromPath(params.trait)
    //input = intergrate_scores.out
      //  .map{btrait, sig, cor, pred, ug_prs_pcs, predicted -> [btrait, ug_prs_pcs]}
        //.combine(trait)
    //input.view()
    //association_analysis(input)
   
}