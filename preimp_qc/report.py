def write_html_report(dirname, basename, qc_tables_list, qc_plots_list):

    # create a complete HTML file
    from datetime import date
    today = date.today()
    date = today.strftime("%B %d, %Y")

    # get plot files
    # pre_man_plt, pre_qq_plt = dirname + 'ManpreQC.png', dirname + 'QQplotpreQC.png'
    # pre_var_cas, pre_var_con = dirname + 'casSNPcrpreQC.png', dirname + 'conSNPcrpreQC.png'
    # pre_id_cas, pre_id_con = dirname + 'casIDcrpreQC.png', dirname + 'conIDPcrpreQC.png'

    abstract = '''This is an automatic output of the QC-Step of the preimp-QC-pipeline, created at MGH, December 2010. It is now in version
    <br/>0.1.0. It is supposed to check and clean a GWAS dataset for technical problems and/or uncontrollable population stratification.
    <br/>This script is still under construction, so please be patient with formatting problems. If you have some ideas, questions,
    <br />requests, feel free to write me.
    '''

    # qc_plots_list = [pre_man_qq_base64, pos_man_qq_base64, pre_con_id_base64, pre_cas_id_base64, pos_con_id_base64, pos_cas_id_base64,
    #                      f_stat_plot, pre_con_var_base64, pre_cas_var_base64, pos_con_var_base64, pos_cas_var_base64]

    # qc_tables_list = [size_of_sample_html, exlusion_overview_html]

    text = '''
    <!doctype html>

    <html lang="en">
    <head style="margin:100px;padding:35px">
      <meta charset="utf-8">

      <title>preimpQC Report</title>
      <meta name="description" content="My description">
      <meta name="Lindo Nkambule" content="My author">

      <style>
        .outer-container {
          display: grid;
          grid-template-columns: 1fr 1fr;
        }

        img {
          max-width: 100%;
        }
      </style>

    </head>

    <body style="margin:100px;padding:35px">
      <h1 style="text-align:center;font-weight:normal">QC-Report of ''' + basename + ''' </h1>
      <h2 style="text-align:center;font-weight:normal">Psychiatric Genomics Consortium</h2>
      <h3 style="text-align:center;font-weight:normal;">''' + date + ''' </h3><br>
        <br>

      <h3 style="text-align:center">Abstract</h3>
      <p style="text-align:center;align=center">''' + abstract + ''' </p>
      <p style="text-align:center;align=center">Lindo Nkambule, lnkambul (at) broadinstitute.org</p>

      <h2>Contents</h2>
      <h3>1 Flags</h3>
      <h3>2 General Info</h3>
      <h3 style="padding-left: 2em;font-weight:normal">2.1 Size of Sample</h3>
      <h3 style="padding-left: 2em;font-weight:normal">2.2 Exclusion overview</h3>
      <h3>3 Manhattan</h3>
      <h3 style="padding-left: 2em;font-weight:normal">3.1 Manhattan-Plot - pre-QC</h3>
      <h3 style="padding-left: 2em;font-weight:normal">3.2 Manhattan-Plot - post-QC</h3>
      <h3>4 Per Individual Characteristic Analysis</h3>
      <h3 style="padding-left: 2em;font-weight:normal">4.1 Missing Rates - pre-QC</h3>
      <h3 style="padding-left: 2em;font-weight:normal">4.2 Missing Rates - post-QC</h3>
      <h3 style="padding-left: 2em;font-weight:normal">4.3 Fhet - pre-QC</h3>
      <h3>5 Per SNP Characteristics Analysis</h3>
      <h3 style="padding-left: 2em;font-weight:normal">5.1 pre-QC missing rate</h3>
      <h3 style="padding-left: 2em;font-weight:normal">5.2 pre-QC missing rate</h3>
      <h3>5 Top-Plots for SNPs</h3>

      <h2>1 Flags</h2>
      <h2>1 There should be a table here. It's coming...</h2>

      <h2>2 General Info</h2>
      <h3>2.1 Size of sample General Info</h3>
      ''' + qc_tables_list[0] + '''
      <h3>2.2 Exclusion overview</h3>
      ''' + qc_tables_list[1] + '''

      <h2>3 Manhattan</h2>
      <h3>3.1 Manhattan-Plot - pre-QC</h2>
      ''' + qc_plots_list[0] + '''
      <h3>3.2 Manhattan-Plot - post-QC</h2>
      ''' + qc_plots_list[1] + '''

      <h2>4. Per Individual Characteristics Analysis</h2>
      <h3>4.1. Missing Rates - pre-QC</h3>
      <div class="outer-container">
        <div>
          ''' + qc_plots_list[2] + '''
        </div>
        <div>
          ''' + qc_plots_list[3] + '''
        </div>
      </div>
      
      <h3>4.2. Missing Rates - post-QC</h3>
      <div class="outer-container">
        <div>
          ''' + qc_plots_list[4] + '''
        </div>
        <div>
          ''' + qc_plots_list[5] + '''
        </div>
      </div>
      
      <h3>4.3. Fstat - pre-QC</h3>
        ''' + qc_plots_list[6] + '''
        
      <h2>5. Per SNP Characteristics Analysis</h3>
      <h3>5.1. pre-QC missing rate</h3>
      <div class="outer-container">
        <div>
          ''' + qc_plots_list[7] + '''
        </div>
        <div>
          ''' + qc_plots_list[8] + '''
        </div>
      </div>
      
      <h3>5.2. post-QC missing rate</h3>
      <div class="outer-container">
        <div>
          ''' + qc_plots_list[9] + '''
        </div>
        <div>
          ''' + qc_plots_list[10] + '''
        </div>
      </div>
      
    </body>
    </html>
    '''

    outhtml = dirname + 'report.html'
    file = open(outhtml, "w")
    file.write(text)
    file.close()