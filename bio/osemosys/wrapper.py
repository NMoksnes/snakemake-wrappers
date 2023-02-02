__author__ = "Nandi Moksnes"
__copyright__ = "Copyright 2023, Nandi Moksnes"
__email__ = "nandi@kth.se"
__license__ = "MIT"

#import snakemake
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

shell("glpsol -m {snakemake.input.model} -d {snakemake.input.data} --wlp {snakemake.output} --check")

# #Pre-processing script MODEX ------ KEEP
# rule modify_model_file:
#     message: "Adding MODEX sets to model file"
#     input:
#         "temp/{scenario}/model_{model_run}.txt"
#     output:
#         temp("temp/{scenario}/model_{model_run}_modex.txt")
#     group: "gen_lp"
#     threads:
#         1
#     conda: "../envs/otoole.yaml"
#     shell:
#         "python workflow/scripts/add_modex_sets.py otoole {input} {output}"



# #Creating LP file ----- KEEP
# rule generate_lp_file:
#     message: "Generating the LP file for '{output}'"
#     input:
#         data="temp/{scenario}/model_{model_run}_modex.txt",
#         model=config['model_file']
#     output:
#         temp(expand("temp/{{scenario}}/model_{{model_run}}.lp{zip_extension}", zip_extension=ZIP))
#     benchmark:
#         "benchmarks/gen_lp/{scenario}_{model_run}.tsv"
#     log:
#         "results/log/glpsol_{scenario}_{model_run}.log"
#     conda: "../envs/osemosys.yaml"
#     group: "gen_lp"
#     threads:
#         1
#     shell:
#         "glpsol -m {input.model} -d {input.data} --wlp {output} --check > {log}"

# rule unzip:
#     message: "Unzipping LP file"
#     input:
#         "temp/{scenario}/model_{model_run}.lp.gz"
#     group:
#         "solve"
#     output:
#         temp("temp/{scenario}/model_{model_run}.lp")
#     shell:
#         "gunzip -fcq {input} > {output}"

# rule solve_lp:
#     message: "Solving the LP for '{output}' using {config[solver]}"
#     input:
#         "temp/{scenario}/model_{model_run}.lp"
#     output:
#         json="modelruns/{scenario}/model_{model_run}/{model_run}.json",
#         solution=temp("temp/{scenario}/model_{model_run}.sol")
#     log:
#         "results/log/solver_{scenario}_{model_run}.log"
#     params:
#         ilp="results/{scenario}/model_{model_run}/solve.ilp",
#         cplex="results/{scenario}/model_{model_run}/solve.cplex",
#     benchmark:
#         "benchmarks/solver/{scenario}_{model_run}.tsv"
#     resources:
#         mem_mb=30000,
#         disk_mb=20000,
#         time=720
#     group: "solve"
#     threads:
#         3
#     shell:
#         """
#         if [ {config[solver]} = gurobi ]
#         then
#           gurobi_cl Method=2 Threads={threads} LogFile={log} LogToConsole=0 ScaleFlag=2 NumericFocus=3 ResultFile={output.solution} ResultFile={output.json} ResultFile={params.ilp} {input}
#         elif [ {config[solver]} = cplex ]
#         then
#           echo "set threads {threads}"   > {params.cplex}
#           echo "set timelimit 43200"     >> {params.cplex}
#           echo "read {input}" 	         >> {params.cplex}
#           echo "baropt"                  >> {params.cplex}
#           echo "write {output.solution}" >> {params.cplex}
#           echo "quit"                    >> {params.cplex}
#         cplex < {params.cplex} > {log} && touch {output.json}
#         else
#           cbc {input} solve -sec 1500 -solu {output.solution} 2> {log} && touch {output.json}
#         fi
#         """

# rule zip_solution:
#     message: "Zip up solution file {input}"
#     group: "solve"
#     input: "temp/{scenario}/model_{model_run}.sol"
#     output: expand("temp/{{scenario}}/{{model_run}}.sol{zip_extension}", zip_extension=ZIP)
#     shell: "gzip -fcq {input} > {output}"

# rule unzip_solution:
#     message: "Unzip solution file {input}"
#     group: "results"
#     input: "temp/{scenario}/model_{model_run}.sol.gz"
#     output: temp("temp/{scenario}/model_{model_run}.sol")
#     shell: "gunzip -fcq {input} > {output}"

# rule transform_file:
#     message: "Transforming CPLEX sol file '{input}'"
#     group: 'results'
#     input: rules.unzip_solution.output
#     conda: "../envs/otoole.yaml"
#     output:
#         temp("temp/{scenario}/model_{model_run}_trans.sol")
#     shell:
#         "python workflow/scripts/transform_31072013.py {input} {output}"

# rule sort_transformed_solution:
#     message: "Sorting transformed CPLEX sol file '{input}'"
#     group: 'results'
#     input:
#         "temp/{scenario}/model_{model_run}_trans.sol"
#     output:
#         temp("temp/{scenario}/model_{model_run}_sorted.sol")
#     shell:
#         "sort {input} > {output}"

