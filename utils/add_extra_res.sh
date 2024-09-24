print_banner(){
    echo "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
}

print_banner
echo "Checking for ligand setting..."
print_banner
# String containing a list of ligand names separated by the pipe symbol "|"
# ligand_fas="FAD|RTC|HEM|DDC|../../mogro/ligand/UDG/UDG.fa.params"
# usage:
# ligand_fas="FAD|RTC|HEM|DDC|../../mogro/ligand/UDG/UDG.fa.params" source $PROTEIN_DESIGN_KIT/4._Tools/scripts/module/add_extra_res.sh

# Initialize an empty string to store the ligand flags
ligand_flags=''

# Split the ligand string into an array
IFS='|' read -ra ligand_array <<< "$ligand_fas"

# Loop through the elements in the array
for i in "${ligand_array[@]}"; do
    # Output the current ligand name
    #
    print_banner
    echo Adding ligand parameters: "$i"
    ligand_pth=""
    # Check if the length of the current ligand name is greater than 3 characters
    if [[ ${#i} -gt 3 && -f ${i} ]];then
      # If so, interpret it as a custom parameter and resolve its absolute path
      ligand_pth=$(readlink -f "$i")
    else
        # If not, search for the ligand in various directories
        for posible_ligand_dir in $(find ./${i}/ ./ligand/${i}/ ./ligands/${i}/ ${PROTEIN_DESIGN_KIT}/5_DB/ligands/${i} -name "${i}.fa.params"); do
          # Check if the file exists
          if [[ -f $posible_ligand_dir ]]; then
            # If it does, output the found location and set it as the path to the ligand
            echo "-->find $i in $posible_ligand_dir"
            ligand_pth=$posible_ligand_dir
            break
          fi
        done
    fi

    # Check if the path to the ligand has been found
    if [[ "$ligand_pth" != "" && -f $ligand_pth ]];then
        # If it has, add the ligand flag to the `ligand_flags` string
        ligand_flags=$ligand_flags"-extra_res_fa $ligand_pth "
    else
        # If the ligand was not found, output an error message
        echo "extra_res not found:  $i $ligand_pth"
    fi
done

print_banner
# Output the final ligand flags string
echo -e "ligand_flags: \n $ligand_flags"

print_banner
