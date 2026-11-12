process generate_xmls {

  input:
    path alignment
    path templates
    val prefix

  output:
    path "*.xml"

  script:
    extra = ""
    if ( params.template_sg )
        extra += " --template-sg ${file(params.template_sg).name}"
    if ( params.template_mat )
        extra += " --template-mat ${file(params.template_mat).name}"
    if ( params.template_egc )
        extra += " --template-egc ${file(params.template_egc).name}"
    if ( params.template_lgc )
        extra += " --template-lgc ${file(params.template_lgc).name}"
    if ( params.template_sg || params.template_mat )
        extra += " --anchor-date ${params.anchor_date}"
    """
    n=\$(grep -c '>' ${alignment})
    generate_beast_xml.py -a ${alignment} -o . -p ${prefix}.n\$n \
        --chain-length ${params.chain_length} --log-every ${params.log_every} ${extra}
    """
}
