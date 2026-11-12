// BDBV phylodynamics pipeline (Part 2): BEAST XML generation, and postprocessing of BEAST output
// (log combining with an ESS gate, tree burnin removal, resampling and MCC tree).
// Nextflow port of bin/run.py. BEAST itself is run outside this pipeline between the two stages.

include { generate_xmls } from './modules/beast'
include { postprocess_logs; postprocess_trees } from './modules/postprocess'

workflow {
  if ( !params.alignment && !params.beast_dir )
      error "Provide --alignment (generate BEAST XMLs) and/or --beast_dir (postprocess BEAST output)"

  // Stage 1: BEAST XMLs from the final alignment of Part 1
  if ( params.alignment ) {
      template_files = [params.template_sg, params.template_mat, params.template_egc, params.template_lgc].findAll { it }
      if ( !template_files )
          error "No templates given: set at least one of --template_sg, --template_mat, --template_egc, --template_lgc"
      alignment_ch = channel.value(file("${params.alignment}", checkIfExists:true))
      templates_ch = channel.fromPath(template_files, checkIfExists:true).collect()
      generate_xmls(alignment_ch, templates_ch, params.prefix)
  }

  // Stage 2: postprocess BEAST output
  if ( params.beast_dir ) {
      logs_ch = channel.fromPath("${params.beast_dir}/*.log", checkIfExists:true)
          .filter { f -> !(f.name =~ /\.(rates|combined)\.log$/) }
          .collect()
      trees_ch = channel.fromPath("${params.beast_dir}/*.trees", checkIfExists:true)
          .filter { f -> !(f.name =~ /burninRemoved|resampled|history|default|empirical/) }
          .collect()

      postprocess_logs(logs_ch)

      // ESS gate: only continue to the trees when ESS is fine, or --force_low_ess is given
      ess_ok = postprocess_logs.out.status
          .filter { f -> params.force_low_ess || f.text.trim() == 'pass' }
          .map { f -> f.text.trim() }
      postprocess_trees(trees_ch, ess_ok)
  }
}

// ESS gate message, printed on completion so it is the last thing on screen.
workflow.onComplete {
  if ( !workflow.success || !params.beast_dir || params.force_low_ess )
      return
  def status = file("${params.outdir}/postprocess_logs/ess_status.txt")
  if ( status.exists() && status.text.trim() == 'low' ) {
      log.info """
      ==========================================================================
      Stopped at the ESS gate: key parameters have ESS below ${params.min_ess}.

      See ${params.outdir}/postprocess_logs/${params.prefix}.ess_report.txt

      Options:
        1. Run more BEAST chains, add their logs and trees to ${params.beast_dir}, and rerun
        2. Continue anyway:

           ${workflow.commandLine} --force_low_ess -resume
      ==========================================================================
      """.stripIndent()
  }
}
