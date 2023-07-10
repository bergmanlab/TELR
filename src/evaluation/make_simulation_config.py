import sys
import os
import json

reference_genomes = sys.argv[1:3]
output = sys.argv[3]

cov_range = [50,100,150,200]

ploidy_dict = {
    2:{
        "name":"diploid",
        "genotypes":{
            "heterozygous":{
                "cov1": lambda cov: cov/2,
                "cov2": lambda cov: cov/2
            },
            "homozygous":{
                "cov1": lambda cov: 0,
                "cov2": lambda cov: cov
            }
        }
    },
    4:{
        "name":"tetraploid",
        "genotypes":{
            "simplex":{
                "cov1": lambda cov: cov*3/4,
                "cov2": lambda cov: cov/4
            },
            "duplex":{
                "cov1": lambda cov: cov/2,
                "cov2": lambda cov: cov/2
            },
            "triplex":{
                "cov1": lambda cov: cov/4,
                "cov2": lambda cov: cov*3/4
            },
            "quadruplex":{
                "cov1": lambda cov: 0,
                "cov2": lambda cov: cov,
            }
        }
    }
}

out_dict = {
    "simulation parameters":{}
}

for ploidy in ploidy_dict:
    ploidy_genotypes = ploidy_dict[ploidy]["genotypes"]
    for genotype in ploidy_genotypes:
        for cov in cov_range:
            dataset = f"{cov}x_{ploidy_dict[ploidy]['name']}_{genotype}"
            out_dict["simulation parameters"][dataset] = {
                os.path.basename(reference_genomes[x]):{
                    "fa":reference_genomes[x],
                    "cov":ploidy_genotypes[genotype][f"cov{x+1}"](cov)
                } for x in [0,1]
            }

with open(output, "w") as out:
    json.dump(out_dict, out)

'''

# run TELR on synthetic tetraploid simplex dataset
for cov in 50 100 150 200; do
    cov1=$((cov * 3 / 4))
    cov2=$((cov / 4))
    ploidy="4"
    genotype="simplex"
    prefix=$cov"x_tetraploid_simplex"
    sample_dir=$out_dir/$prefix
    log_output=$out_dir/$prefix.log
    sbatch --job-name=${prefix} --export=utility_dir="${utility_dir}",region1_mask="${a4_region_mask}",out_dir="${sample_dir}",reads="${fq}",ref1="${fa1}",ref2="${fa2}",cov1="${cov1}",cov2="${cov2}",prefix="${prefix}",cov2="${cov2}",te_library="${te_library}",ref1_te_annotation="${lift_nonref_te}",ref2_te_annotation="${dm6_te}",ploidy="${ploidy}",genotype="${genotype}" --output="${log_output}" --error="${log_output}" $eval_slurm_job
done

# run TELR on synthetic tetraploid duplex dataset
for cov in 50 100 150 200; do
    cov1=$((cov / 2))
    cov2=$((cov / 2))
    ploidy="4"
    genotype="duplex"
    prefix=$cov"x_tetraploid_duplex"
    sample_dir=$out_dir/$prefix
    log_output=$out_dir/$prefix.log
    sbatch --job-name=${prefix} --export=utility_dir="${utility_dir}",region1_mask="${a4_region_mask}",out_dir="${sample_dir}",reads="${fq}",ref1="${fa1}",ref2="${fa2}",cov1="${cov1}",cov2="${cov2}",prefix="${prefix}",cov2="${cov2}",te_library="${te_library}",ref1_te_annotation="${lift_nonref_te}",ref2_te_annotation="${dm6_te}",ploidy="${ploidy}",genotype="${genotype}" --output="${log_output}" --error="${log_output}" $eval_slurm_job
done

# run TELR on synthetic tetraploid triplex dataset
for cov in 50 100 150 200; do
    cov1=$((cov / 4))
    cov2=$((cov * 3 / 4))
    ploidy="4"
    genotype="triplex"
    prefix=$cov"x_tetraploid_triplex"
    sample_dir=$out_dir/$prefix
    log_output=$out_dir/$prefix.log
    sbatch --job-name=${prefix} --export=utility_dir="${utility_dir}",region1_mask="${a4_region_mask}",out_dir="${sample_dir}",reads="${fq}",ref1="${fa1}",ref2="${fa2}",cov1="${cov1}",cov2="${cov2}",prefix="${prefix}",cov2="${cov2}",te_library="${te_library}",ref1_te_annotation="${lift_nonref_te}",ref2_te_annotation="${dm6_te}",ploidy="${ploidy}",genotype="${genotype}" --output="${log_output}" --error="${log_output}" $eval_slurm_job
done

# run TELR on synthetic tetraploid quadruplex dataset
for cov in 50 100 150 200; do
    cov1=0
    cov2=$cov
    ploidy="4"
    genotype="quadruplex"
    prefix=$cov"x_tetraploid_quadruplex"
    sample_dir=$out_dir/$prefix
    log_output=$out_dir/$prefix.log
    sbatch --job-name=${prefix} --export=utility_dir="${utility_dir}",region1_mask="${a4_region_mask}",out_dir="${sample_dir}",reads="${fq}",ref1="${fa1}",ref2="${fa2}",cov1="${cov1}",cov2="${cov2}",prefix="${prefix}",cov2="${cov2}",te_library="${te_library}",ref1_te_annotation="${lift_nonref_te}",ref2_te_annotation="${dm6_te}",ploidy="${ploidy}",genotype="${genotype}" --output="${log_output}" --error="${log_output}" $eval_slurm_job
done

# run TELR on synthetic diploid heterozyguous dataset
for cov in 50 100 150 200; do
    cov1=$((cov / 2))
    cov2=$((cov / 2))
    ploidy="2"
    genotype="heterozygous"

    prefix=$cov"x_diploid_heterozygous"
    sample_dir=$out_dir/$prefix
    log_output=$out_dir/$prefix.log
    sbatch --job-name=${prefix} --export=utility_dir="${utility_dir}",region1_mask="${a4_region_mask}",out_dir="${sample_dir}",reads="${fq}",ref1="${fa1}",ref2="${fa2}",cov1="${cov1}",cov2="${cov2}",prefix="${prefix}",cov2="${cov2}",te_library="${te_library}",ref1_te_annotation="${lift_nonref_te}",ref2_te_annotation="${dm6_te}",ploidy="${ploidy}",genotype="${genotype}" --output="${log_output}" --error="${log_output}" $eval_slurm_job
done

# run TELR on synthetic diploid homozygous dataset
for cov in 50 100 150 200; do
    cov1=0
    cov2=$cov
    ploidy="2"
    genotype="homozygous"

    prefix=$cov"x_diploid_homozygous"
    sample_dir=$out_dir/$prefix
    log_output=$out_dir/$prefix.log
    sbatch --job-name=${prefix} --export=utility_dir="${utility_dir}",region1_mask="${a4_region_mask}",out_dir="${sample_dir}",reads="${fq}",ref1="${fa1}",ref2="${fa2}",cov1="${cov1}",cov2="${cov2}",prefix="${prefix}",cov2="${cov2}",te_library="${te_library}",ref1_te_annotation="${lift_nonref_te}",ref2_te_annotation="${dm6_te}",ploidy="${ploidy}",genotype="${genotype}" --output="${log_output}" --error="${log_output}" $eval_slurm_job
done
'''