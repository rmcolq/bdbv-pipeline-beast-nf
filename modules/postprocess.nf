// Combine chain logs, remove burnin and check ESS. A low ESS does not fail the process:
// it is recorded in ess_status.txt (pass / low) so the workflow can stop at the ESS gate.
process postprocess_logs {

  input:
    path logs

  output:
    path "${params.prefix}.combined.log", emit: combined
    path "${params.prefix}.ess_report.txt", emit: report
    path "ess_status.txt", emit: status

  script:
    """
    rc=0
    postprocess_logs.py --logs ${logs} --burnin ${params.burnin} --outdir . \
        --prefix ${params.prefix} --min-ess ${params.min_ess} || rc=\$?
    if [ \$rc -eq 0 ]; then
        echo pass > ess_status.txt
    elif [ \$rc -eq 2 ]; then
        echo low > ess_status.txt
    else
        exit \$rc
    fi
    """
}

process postprocess_trees {

  input:
    path trees
    val ess_ok

  output:
    path "${params.prefix}.burninRemoved.trees"
    path "${params.prefix}.burninRemoved.resampled.trees"
    path "${params.prefix}.mcc.tree", optional: true

  script:
    """
    postprocess_trees.py --trees ${trees} --burnin ${params.burnin} --resample ${params.resample} \
        --outdir . --prefix ${params.prefix} --treeannotator ${params.treeannotator}
    """
}
