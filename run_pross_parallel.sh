#! /bin/bash
set -e
## arranged by Yinying Yao for Rosetta ver. 3.12
# Scorefunction talaris 2014 is replaced with ref2015

# PROSS
# Publication describing the method:
# Adi Goldenzweig, Moshe Goldsmith, Shannon E Hill, Or Gertman, Paola Laurino, Yacov Ashani, Orly Dym, Tamar Unger, Shira Albeck, Jaime Prilusky, Raquel L Lieberman, Amir Aharoni, Israel Silman, Joel L Sussman, Dan S Tawfik and Sarel J Fleishman
# Automated Structure- and Sequence-Based Design of Proteins for High Bacterial Expression and Stability
# Molecular cell (published online on July 14 2016)

script_dir=$(dirname $0)
MODULE_PATH=${script_dir}/utils/

if ! command -v  whichrosetta; then
    pip install rosetta_finder -U
fi

if ! command -v  parallel; then
    echo 'GNU-parallel is required for run this PROSS protocol.'
    exit 1
fi

# USAGE

# # FETCH OPTIONS IF NOT SLURM
# ####################################################################################################################
# if [[ $IN_SLURM == false ]]; then
usage() {
    echo ""
    echo "Usage: $0 <OPTIONS>"
    echo "Required Parameters:"
    echo "      -s                  <init_pdb> "
    echo "Optional Parameters:"
    echo "      -f                  <res_to_fix> "
    echo "      -r                  <res_to_restrict> "
    echo "      -p                  <pssm_fn> "
    echo "      -l                  <ligand_fas> "
    echo "      -j                  <nproc> "
    echo "      -c                  <chain_id> "
    echo "      -d   true|false     <delay> "
    echo "      -o                  <output_dir> "
    echo "      -u                  <uniref90_db> Path/prefix to Uniref90, default for JAPS sever. "
    echo "      -B                  <blast_bin> Path/prefix to NCBI BLAST, default as $(dirname $(which psiblast))"
    echo "      -h                  print help message and exit."
    echo ""
    exit 1
}
while getopts ":s:p:j:l:f:r:c:o:u:B:dh" opt; do
    case "${opt}" in
    # required options
    s) init_pdb=$OPTARG ;;
    p) pssm_fn=$OPTARG ;;
    f) res_to_fix=$OPTARG ;;
    r) res_to_restrict=$OPTARG ;;
    l) ligand_fas=$OPTARG ;;
    j) nproc=$OPTARG ;;
    c) chain_id=$OPTARG ;;
    d) delay=true ;;
    o) output_dir=$OPTARG ;;
    u) uniref90_db=$OPTARG ;;
    B) blast_bin=$OPTARG ;;
    h) usage ;;
    *)
        echo Unknown option!
        usage
        ;;
    esac
done
# preset configs
program_path=/software/blast-2.2.27/ncbi-blast-2.2.27+/bin/
database=/mnt/db/uniref90/uniref90
BIN_ROSETTA_SCRIPTS=$(whichrosetta rosetta_scripts)
# else
#     # PSSM config , new slurm cluster
#     if [[ $(hostname) =~ "arm" ]];then
#         program_path=/public/software/env01-arm64/bin/
#     else
#         program_path=/public/agis/huangsanwen_group/yaoyinying/software/ncbi-blast/ncbi-blast-2.2.27+/bin/
#     fi
#     database=/public/agis/huangsanwen_group/yaoyinying/db/uniref90/uniref90
# fi

print_banner(){
    echo "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
}

# FUNCTIONS
SeqLenPDB() {
    local pdb_fp=$1
    local chain_id_referred=$2
    if [[ "$chain_id_referred" == "" ]]; then
        local chain_id_referred=A
    fi
    #local chain_num=$(cat ${pdb_fp} | grep -e "^ATOM" | cut -b 22 |uniq)
    for chain_id in $(cat ${pdb_fp} | grep -e "^ATOM" | cut -b 22 | uniq); do
        if [[ "$chain_id" == "$chain_id_referred" ]]; then
            local chain_len=$(cat ${pdb_fp} | grep -e "^ATOM" | cut -b 22-26 | uniq | grep -e "^${chain_id}" | cut -b 3- | tail -n1)
            echo $chain_len
        fi
    done
    #	return $C
}


RUN_PSSM() {
    local fasta=$(readlink -f $1)
    local fasta_fn=$(basename $1)
    if [[ ! -f $fasta ]]; then
        echo FASTA not exist: $fasta
        exit 1
    fi

    local instance=${fasta_fn%.fasta}

    echo Running BLAST sequence searching ...
    echo Processing $fasta ...

    # expected output_dir
    local dst=$(dirname "${fasta}")
    # expected msa output_dir
    local dst_msa=$dst/pssm_msa

    mkdir -p $dst_msa

    #pushd $dst_msa
    if [[ ! -f "${instance}_ascii_mtx_file" ]]; then
        local cmd="${blast_bin}/psiblast -query ${fasta} -db ${database} -out_pssm ${dst_msa}/${instance}.ckp \
                -evalue 0.001 -out_ascii_pssm ${dst_msa}/${instance}_ascii_mtx_file -out ${dst_msa}/${instance}_output_file \
                -num_iterations 3 -num_threads ${nproc}"
        echo "$cmd"
        eval "$cmd"
    fi
    wait
    #popd

}

# DEFINE SUBDIR
# if the structure of directory changes someday, we need to update the statements here.
path_to_script=$script_dir
echo script path: "${path_to_script}"

# set default parameters
if [[ "$nproc" == "" ]]; then
    nproc=$(nproc)
fi
# PROCESSING UNDEFINED ARGUMENTS
echo Using $nproc processors.

# set default parameters
if [[ "$chain_id" == "" ]]; then
    chain_id=A
fi

if [[ "$delay" == "" ]]; then
    delay=false
fi

if [[ "$delay" == true ]]; then
    flags_pross=flags_delay
else
    flags_pross=flags_nodelay
fi


if [[ ! -f "$init_pdb" ]]; then
    echo required input missing!
    usage
fi

if [[  "$res_to_fix" == ""  ]]; then
    res_to_fix=1A
fi

if [[  "$res_to_restrict" == ""  ]]; then
    res_to_restrict=1A
fi

if [[ "$(python -c "import Bio" 2>&1)" =~ "ModuleNotFoundError" ]]; then
    echo biopython not available. exiting ...
    exit 1
fi

# if given $blast_bin is blank or unavailable
if [[ "$blast_bin" == "" || ! -f "$blast_bin" ]];then
  # try to use default NCBI BLAST+ if its found in PATH
  if [[ "$(which psiblast)" != "" ]];then
    blast_bin=$(dirname $(which psiblast))
  else
    # otherwise, use hardcoded program path
    blast_bin=$program_path
  fi
fi

# if given $uniref90_db is blank or unavailable
if [[ "$uniref90_db" == "" || ! -f "$uniref90_db" ]]; then
    uniref90_db=$database
fi

# SET JOB NAME
if [[ "$job_name" == "" ]]; then
    job_name="PROSS"
fi
echo Using Job Name $job_name.

# PROCESS LIGAND PARAMETERS
source "$MODULE_PATH/add_extra_res.sh"
echo Ligand Parameter: $ligand_fas

# PROCESS ADVANCED ROSETTA FLAGS
source $MODULE_PATH/add_adv_flags.sh

echo Initial model: $init_pdb

pdb_filename=$(basename $init_pdb)
instance=${pdb_filename%.pdb}
initial_model=$(readlink -f "$init_pdb")

if [[ $output_dir == "" ]];then
  output_dir=${job_name}_${instance}_output
fi

mkdir -p ${output_dir}
mkdir -p ${output_dir}/logs

print_banner
echo "Step 0: prepare PSSM profile."
if [[ "$pssm_fn" == "" ]]; then
    if [[ ! -f ${output_dir}/pssm_msa/${instance}_ascii_mtx_file ]]; then
        fasta_fn=${instance}.fasta
        echo '>'${instance} >${output_dir}/${instance}.fasta
        python3 ${MODULE_PATH}/sequence_extract.py ${initial_model} A >>${output_dir}/${instance}.fasta

        # expected sequence
        RUN_PSSM ${output_dir}/${instance}.fasta
        # expected msa: ${output_dir}/pssm_msa/${instance}_ascii_mtx_file
    else
        echo PSSM file found: ${output_dir}/pssm_msa/${instance}_ascii_mtx_file
    fi

    # expected pssm output: $(output_dir)/$(instance)/pssm_msa/${instance}_ascii_mtx_file
    if [[ ! -f ${output_dir}/pssm_msa/${instance}_ascii_mtx_file ]]; then
        echo PSSM file not found: ${output_dir}/pssm_msa/${instance}_ascii_mtx_file
        echo RUN_PSSM may be failed.
        exit 1
    else
        pssm_fn=${output_dir}/pssm_msa/${instance}_ascii_mtx_file
    fi
else
    pssm_fn=$(readlink -f $pssm_fn)
fi


# BUILD COMMAND AND RUN IT
print_banner
echo "Step 1: Build the CST file from init_pdb."
# step 1:
# make contraints
cst_file="${instance}_bbCA.cst"
if [[ ! -f ${output_dir}/${cst_file}  || $(cat ${output_dir}/${cst_file} |wc -l |awk '{print $NF}') != $(SeqLenPDB ${init_pdb} A ) ]];then
    echo CST file not found or uncomplete.
    cmd="awk -f  ${path_to_script}/scripts/make_cst.awk ${initial_model} > ${output_dir}/${cst_file}"
    echo "$cmd"
    eval "$cmd"
else
    echo cst_file found: ${output_dir}/${cst_file}
fi

# step 2:
# structure refinement
print_banner
echo "Step 2: Structure refinement"

refine_res_dir=${output_dir}/refinement/
mkdir -p ${refine_res_dir}
if [[ ! -f ${refine_res_dir}/${instance}.pdb ]];then
    cmd=" ${BIN_ROSETTA_SCRIPTS} \
    -in:file:s ${initial_model} ${ligand_flags} \
    -out:path:all ${refine_res_dir}/ \
    -parser:protocol ${path_to_script}/scripts/refine_new.xml \
    -parser:script_vars res_to_fix=${res_to_fix} \
    -parser:script_vars pdb_reference=${initial_model} \
    -parser:script_vars cst_full_path=${output_dir}/${cst_file}  \
    -parser:script_vars cst_value=0.4 @${path_to_script}/scripts/flags_nodelay \
    -no_nstruct_label true "

    echo "$cmd"
    eval "$cmd" >${output_dir}/logs/${job_name}_${instance}_pross_refine.log
else
    echo "find previous run result of refinement: ${refine_res_dir}/${instance}.pdb"
fi

print_banner
echo "Step 3: Run the filterscan"

filterscan_res_dir=${output_dir}/filterscan/
mkdir -p ${filterscan_res_dir}
mkdir -p ${filterscan_res_dir}/resfiles
mkdir -p ${filterscan_res_dir}/scores
mkdir -p ${filterscan_res_dir}/resfiles/tmp
if [[ -f "${filterscan_res_dir}/resfiles/designable_aa_resfile.0.5" || -f "${filterscan_res_dir}/resfiles/designable_aa_resfile.-0.45" || -f "${filterscan_res_dir}/resfiles/designable_aa_resfile.-1" ]];then
    echo Previous filterscan file fond. Now skipping this step.
else
    seq_length=$(SeqLenPDB ${init_pdb} ${chain_id})
    
    parallel --results ${output_dir}/logs/filterscan_parallel/${instance}/ --jobs ${nproc} -k --bar -v \
    echo Running filterscan for residue: {} ';' \
    ${BIN_ROSETTA_SCRIPTS} \
    -in:file:s ${refine_res_dir}/${instance}.pdb ${ligand_flags} \
    -out:path:all ${filterscan_res_dir}  \
    -out:prefix ${instance}_pross_filterscan_ \
    -out:file:scorefile ${instance}_pross_filterscan.sc \
    -parser:protocol ${path_to_script}/scripts/filterscan_parallel.xml \
    -parser:script_vars res_to_fix=${res_to_fix} \
    -parser:script_vars pdb_reference=${refine_res_dir}/${instance}.pdb \
    -parser:script_vars res_to_restrict=${res_to_restrict} \
    -parser:script_vars cst_full_path=${output_dir}/${cst_file} \
    -parser:script_vars cst_value=0.4 \
    -parser:script_vars pssm_full_path=${pssm_fn} \
    -parser:script_vars scores_path=${filterscan_res_dir}/scores/ \
    -parser:script_vars resfiles_path=${filterscan_res_dir}/resfiles/tmp/ \
    @${path_to_script}/scripts/${flags_pross} \
    -parser:script_vars current_res={} -overwrite ::: $(seq 1 ${seq_length})

    

    # expected resfile: ${filterscan_res_dir}/resfiles/tmp/designable_aa_resfile-${res_id}.${level}

    print_banner
    echo "Step 4 :merge resfiles."
    for level in 0.5 -0.45 -0.75 -1 -1.25 -1.5 -1.8 -2;do
        resfile_fn=designable_aa_resfile.${level}
        first_resfile=true
        for res_id in $(seq 1 ${seq_length});do
            tmp_resfile_fn=designable_aa_resfile-${res_id}.${level}
            if [[ ! -f ${filterscan_res_dir}/resfiles/tmp/${tmp_resfile_fn} ]];then
                echo TmpResFile not found: ${tmp_resfile_fn}
                continue;
            elif [[ $first_resfile == true ]];then
                echo  Head TmpResFile found: ${tmp_resfile_fn}
                cat ${filterscan_res_dir}/resfiles/tmp/${tmp_resfile_fn} >${filterscan_res_dir}/resfiles/${resfile_fn}
                first_resfile=false
            else
                # append to the resfile
                cat ${filterscan_res_dir}/resfiles/tmp/${tmp_resfile_fn} |grep -e "^[[:digit:]]" >>${filterscan_res_dir}/resfiles/${resfile_fn}
            fi
        done
    done
fi
# 4. design


print_banner
echo "Step 5: make designs."
design_res_dir=${output_dir}/design/
mkdir -p ${design_res_dir}
mkdir -p ${design_res_dir}/scores


resfile_list=$(ls ${filterscan_res_dir}/resfiles/ |grep designable_aa_resfile)
echo "Run design to these resfiles: "
echo ${resfile_list}
parallel --plus  --results ${output_dir}/logs/design/${instance}/ --jobs ${nproc} -k --bar -v  \
    ${BIN_ROSETTA_SCRIPTS} \
    -in:file:s ${refine_res_dir}/${instance}.pdb ${ligand_flags} \
    -out:path:pdb ${design_res_dir} \
    -out:path:score ${design_res_dir}/scores \
    -out:suffix _{#designable_aa_resfile.} \
    -no_nstruct_label true \
    -out:file:scorefile ${instance}_pross_design_{#designable_aa_resfile.}.sc \
    -parser:protocol ${path_to_script}/scripts/design_new.xml \
    -parser:script_vars in_resfile=${filterscan_res_dir}/resfiles/{} \
    -parser:script_vars res_to_fix=${res_to_fix} \
    -parser:script_vars pdb_reference=${refine_res_dir}/${instance}.pdb \
    -parser:script_vars cst_full_path=${output_dir}/${cst_file} \
    -parser:script_vars cst_value=0.4 \
    -parser:script_vars pssm_full_path=${pssm_fn} @${path_to_script}/scripts/flags_nodelay ::: ${resfile_list}

print_banner
echo "PROSS is done! ${instance}"

mkdir -p ${output_dir}/report/

echo "Now generate design report ."
for pdb in ${design_res_dir}/*.pdb;do
  echo $(readlink -f $pdb) | awk '{
      # parse basename of pdb
      pdb_bn_idx=split($1,bn_arr,"/");
      pdb_bn=bn_arr[pdb_bn_idx];

      # parse design level of PROSS coded in pdb filename
      level_idx=split(pdb_bn,level_arr,"_");
      level=level_arr[level_idx];
      sub(".pdb","",level);

      # print into `design_level,absolute-path-of-designer`
      print "PROSS_"level","$1
    }';
done > ${output_dir}/report/${instance}_${job_name}_report.txt

print_banner

