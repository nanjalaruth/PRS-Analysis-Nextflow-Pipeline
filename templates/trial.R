"""
        # Check if the columns chr_name and chr_position exist in the input file
        if awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++){if(\$i=="chr_name"){c1=1}if(\$i=="chr_position"){c2=1}}}END{exit !(c1 && c2)}'\
         ${score_files}; then \
             awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++) { if (\$i == "chr_name") \
              { chr_name_index=i } else if (\$i == "chr_position") \
                { chr_position_index=i } } print \$0, "\\tSNP"; next} \
                 {print \$0 "\\t" \$chr_name_index ":" \$chr_position_index}' \
                   ${score_files} > ${blood_trait}_${pgs_id}_edit_2.txt
        else
            # Replace rsID with chr_name:chr_position using the reference file
            awk -F'\\t' '
            BEGIN {
                while ((getline < "snp_info_header.txt") > 0) {
                    map[\$1] = \$1":"\$2
                }
                close("snp_info_header.txt")
            }
            NR==1 {print \$0; next} {if(\$1 in map) \$1=map[\$1]; print}' ${score_files} > ${blood_trait}_${pgs_id}_edit_2.txt
        fi
        """


        awk -F'\\t' '
            BEGIN {
                while ((getline < "snp_info_header.txt") > 0) {
                    map[\$1] = \$1":"\$2
                }
                close("snp_info_header.txt")
            }
            NR==1 {print \$0; next} {if(\$1 in map) \$1=map[\$1]; print}' ${score_files} > ${blood_trait}_${pgs_id}_edit_2.txt
        fi