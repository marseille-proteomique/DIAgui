library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(shinyWidgets)
library(shinyFiles)
library(data.table)
library(dplyr)
library(plotly)
library(stringr)
library(Rcpp)
library(RcppEigen)
library(iq)


#increase the max request size for uploading files
options(shiny.maxRequestSize = 10000*1024^2)
#set options for the spinner when things are loading
options(spinner.color = "#3D7AEC", spinner.color.background = "000000", spinner.size = 2)

ui <- fluidPage(
  useShinydashboard(), #allow to use box without dashboard

  tags$head(tags$style(HTML(".navbar-default {background-color: #3D7AEC !important; color = #ffffff}
                             .navbar-default > .container-fluid > .navbar-nav > li > a {color: #ffffff; font-size: 18px}
                             .navbar-default > .container-fluid > .navbar-nav > li > a:hover {background-color: #165CDE; color: #ffffff}
                             .navbar-default > .container-fluid > .navbar-nav > li[class=active] > a {background-color: #003EFF; color: #ffffff}
                             .navbar-default > .container-fluid > .navbar-nav > li[class=active] > a:hover {background-color: #003EFF; color: #ffffff}
                             .navbar-default > .container-fluid > .navbar-header > .navbar-brand {color: #ffffff; font-size: 22px}
                             * {font-family: 'Rockwell'}
                             body {background-color: #C0CEFF}
                             .nav-tabs > li > a {background-color: #3D7AEC; color: #ffffff}
                             .nav-tabs > li > a:hover {background-color: #165CDE; color: #ffffff}
                             .nav-tabs > li[class=active] > a {background-color: #003EFF; color: #ffffff}
                             .nav-tabs > li[class=active] > a:hover {background-color: #003EFF; color: #ffffff}
                             .datatables {background-color: #ffffff")
                       )
            ),

  navbarPage(
    title = "DIA-NN R routine",
    tabPanel("Process",
             sidebarLayout(
               sidebarPanel(
                 conditionalPanel(condition = "input.tab1 != 'Import your data' & output.reportdata_up",
                                  HTML("<p><h3>General info</h3><br>
                                           All quantities are based on the column named 'Precursor.Normalized' from the report file.<br>
                                           Each threshold you can select correspond to a q-value (look in the report your imported).
                                           If you set a value to 1, it will not apply any filter according to this value.<br>
                                           All generated files can be saved in txt, csv or xlsx format.
                                        <p><h3>The MaxLFQ algorithm</h3><br>
                                               This algorithm is another way to determine intensity and to normalize data
                                               in Label-Free quantification.  Quickly, the aim is to perfom a 'delayed normalization'
                                               by determining normalization coefficients for each fraction, and then
                                               extracts the maximum ratio information from peptide signals in arbitrary
                                               numbers of samples to achieve the highest possible accuracy of quantification.<br>
                                               For more information, see this <a href=https://pubmed.ncbi.nlm.nih.gov/24942700/>article</a>
                                               from Jurgen Cox and al.
                                        </p>"
                                       )
                                  ),
                 conditionalPanel(condition = "input.tab1 == 'Import your data'",
                                  HTML("<p>After your analysis with the DIA-nn software, you can find the file 'report.tsv' in your results.
                                           From this file, you can filter your data according to your criterias, use the MaxLFQ algorithm
                                           for quantification and normalization, get the number of peptides used for the quantification,
                                           get the Top3 absolute quantification and then save the files you want to keep. <br>
                                           The aim of this app is to provide you a 'user-friendly' interface in order to use the diann R routine
                                           and adding some other usefull informations.
                                        </p>"
                                       )
                                  ),
                 width = 3
                 ),
               mainPanel(
                 tabsetPanel(type = "tabs", id = "tab1",
                             tabPanel("Import your data",
                                      tags$hr(),
                                      fluidRow(box(title = "DIA data", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                   fluidRow(column(6, shinyFilesButton("rep_tsv", label = "Select your report file", title = "Please select a file",
                                                                                       icon = icon("file"),
                                                                                       multiple = TRUE, viewtype = "detail", buttonType = "primary", class = "btn-lg"))
                                                            ),
                                                   tags$hr(),
                                                   conditionalPanel(condition = "output.reportdata_up",
                                                                    tags$u(h2("Your data")),
                                                                    fluidRow(column(12, DT::dataTableOutput("df_report"))
                                                                             ),
                                                                    tags$hr(),
                                                                    tags$u(h2("Rename your fractions")),
                                                                    htmlOutput("frac_dat"),
                                                                    tags$hr(),
                                                                    radioButtons("chorename_dat", "Choose how to rename your fractions",
                                                                                 choices = c("Remove path", "New names"), selected = "Remove path", inline = TRUE),
                                                                    fluidRow(conditionalPanel(condition = "input.chorename_dat == 'New names'",
                                                                                              column(6, textInput("newfrac_dat", "Type the new name of your fraction,
                                                                                                  separated with a comma, in the same order.
                                                                                                  (if empty between comma, no changes)"))
                                                                                              ),
                                                                             conditionalPanel(condition = "input.chorename_dat == 'Remove path'",
                                                                                              column(6, radioButtons("whattorm_dat", "Choose what to remove",
                                                                                                                     choices = list("Only keep what's different" = 1,
                                                                                                                                    "Only keep file name with extension" = 2,
                                                                                                                                    "Only keep file name without extension" = 3),
                                                                                                                     selected = 3, inline = TRUE
                                                                                                                     )
                                                                                                     )
                                                                                              ),
                                                                             column(4, actionButton("change_dat", "Rename your fraction", class = "btn-primary"))
                                                                             )
                                                                    ),
                                                   tags$hr()
                                                   )
                                               )

                                      ),
                             tabPanel("Peptides and precursors",
                                      conditionalPanel(condition = "!output.reportdata_up",
                                                       h2("Import a report file in the tab 'Import your data' first !")
                                                       ),
                                      conditionalPanel(condition = "output.reportdata_up",
                                                       tags$hr(),
                                                       fluidRow(box(title = "Precursors", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                                    tags$u(h3("Get your precursor file")),
                                                                    tags$hr(),

                                                                    fluidRow(column(3, numericInput("qv_prec", "Choose the q-value to filter the precursors",
                                                                                                    min = 0, max = 1, step = 0.01, value = 0.01)),
                                                                             column(3, numericInput("qvprot_prec", "Choose the protein.q-value to filter the precursors",
                                                                                                    min = 0, max = 1, step = 0.01, value = 1)),
                                                                             column(3, numericInput("qvpg_prec", "Choose the protein-group q-value to filter the precursors",
                                                                                                    min = 0, max = 1, step = 0.01, value = 0.01)),
                                                                             column(3, numericInput("qvgg_prec", "Choose the gene-groupe q-value to filter the precursors",
                                                                                                    min = 0, max = 1, step = 0.01, value = 1))
                                                                             ),
                                                                    fluidRow(column(3, checkboxInput("protypiconly_prec", "Proteotypic only", TRUE)),
                                                                             column(3, radioButtons("remmodif_prec", "",
                                                                                                    choices = c("Only keep modification selected" = "keep",
                                                                                                                "Remove modifications selected" = "rem")
                                                                                                    )
                                                                                    ),
                                                                             column(3, selectInput("modif_prec", "Choose some modifications (if NULL, no filtering)",
                                                                                                   choices = "",
                                                                                                   multiple = TRUE))
                                                                             ),
                                                                    actionButton("go_prec", "Start calculation", class = "btn-primary"),
                                                                    tags$hr(),
                                                                    textOutput("info_prec"),
                                                                    conditionalPanel(condition = "output.precursor_up",
                                                                                     DT::dataTableOutput("res_prec"),
                                                                                     tags$hr(),
                                                                                     fluidRow(column(3, downloadButton("down_prec", "Download results")),
                                                                                              column(3, selectInput("format_prec", "Select a format",
                                                                                                                    choices = c("txt", "csv", "xlsx"),
                                                                                                                    selected = "txt"))
                                                                                              )
                                                                                     )
                                                                    )
                                                                ),
                                                       fluidRow(box(title = "Peptides", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                                    tags$u(h3("Get your peptide file")),
                                                                    tags$hr(),

                                                                    fluidRow(column(3, numericInput("qv_pep", "Choose the q-value to filter the peptides",
                                                                                                    min = 0, max = 1, step = 0.01, value = 0.01)),
                                                                             column(3, numericInput("qvprot_pep", "Choose the protein.q-value to filter the peptides",
                                                                                                    min = 0, max = 1, step = 0.01, value = 1)),
                                                                             column(3, numericInput("qvpg_pep", "Choose the protein-group q-value to filter the peptides",
                                                                                                    min = 0, max = 1, step = 0.01, value = 0.01)),
                                                                             column(3, numericInput("qvgg_pep", "Choose the gene-groupe q-value to filter the peptides",
                                                                                                    min = 0, max = 1, step = 0.01, value = 1))
                                                                             ),
                                                                    fluidRow(column(3, checkboxInput("protypiconly_pep", "Proteotypic only", TRUE)),
                                                                             column(3, selectInput("centercol_pep", "Choose identifier to quantify",
                                                                                                   choices = c("Modified.Sequence", "Stripped.Sequence"),
                                                                                                   selected = "Stripped.Sequence")),
                                                                             column(3, radioButtons("remmodif_pep", "",
                                                                                                    choices = c("Only keep modification selected" = "keep",
                                                                                                                "Remove modifications selected" = "rem")
                                                                                                    )
                                                                                    ),
                                                                             column(3, selectInput("modif_pep", "Choose some modifications (if NULL, no filtering)",
                                                                                                   choices = "",
                                                                                                   multiple = TRUE))
                                                                             ),
                                                                    fluidRow(column(3, checkboxInput("getPTM_pep", "Extract best PTM.Q.value", FALSE))
                                                                             ),
                                                                    actionButton("go_pep", "Start calculation", class = "btn-primary"),
                                                                    tags$hr(),
                                                                    textOutput("info_pep"),
                                                                    conditionalPanel(condition = "output.peptide_up",
                                                                                     DT::dataTableOutput("res_pep"),
                                                                                     tags$hr(),
                                                                                     fluidRow(column(3, downloadButton("down_pep", "Download results")),
                                                                                              column(3, selectInput("format_pep", "Select a format",
                                                                                                                    choices = c("txt", "csv", "xlsx"),
                                                                                                                    selected = "txt"))
                                                                                              )
                                                                                     ),

                                                                    tags$u(h3("Get your peptide file using the MaxLFQ algorithm")),
                                                                    tags$hr(),

                                                                    fluidRow(column(3, numericInput("qv_peplfq", "Choose the q-value to filter the peptides",
                                                                                                    min = 0, max = 1, step = 0.01, value = 0.01)),
                                                                             column(3, numericInput("qvprot_peplfq", "Choose the protein.q-value to filter the peptides",
                                                                                                    min = 0, max = 1, step = 0.01, value = 1)),
                                                                             column(3, numericInput("qvpg_peplfq", "Choose the protein-group q-value to filter the peptides",
                                                                                                    min = 0, max = 1, step = 0.01, value = 0.01)),
                                                                             column(3, numericInput("qvgg_peplfq", "Choose the gene-groupe q-value to filter the peptides",
                                                                                                    min = 0, max = 1, step = 0.01, value = 1))
                                                                             ),
                                                                    radioButtons("wLFQ_peplfq", "",
                                                                                 choices = c("Use fast MaxLFQ from iq package (log2 transformed)" = "iq",
                                                                                             "Use MaxLFQ from diann package" = "diann"),
                                                                                 selected = "iq",
                                                                                 inline = TRUE),
                                                                    fluidRow(column(3, checkboxInput("protypiconly_peplfq", "Proteotypic only", TRUE)),
                                                                             column(3, selectInput("centercol_peplfq", "Choose identifier to quantify",
                                                                                                   choices = c("Modified.Sequence", "Stripped.Sequence"),
                                                                                                   selected = "Modified.Sequence")),
                                                                             column(3, radioButtons("remmodif_peplfq", "",
                                                                                                    choices = c("Only keep modification selected" = "keep",
                                                                                                                "Remove modifications selected" = "rem")
                                                                                                    )
                                                                                    ),
                                                                             column(3, selectInput("modif_peplfq", "Choose some modifications (if NULL, no filtering)",
                                                                                                   choices = "",
                                                                                                   multiple = TRUE))
                                                                             ),
                                                                    fluidRow(column(3, checkboxInput("getPTM_peplfq", "Extract best PTM.Q.value", FALSE))
                                                                             ),
                                                                    actionButton("go_peplfq", "Start calculation", class = "btn-primary"),
                                                                    tags$hr(),
                                                                    textOutput("info_peplfq"),
                                                                    conditionalPanel(condition = "output.peptideLFQ_up",
                                                                                     DT::dataTableOutput("res_peplfq"),
                                                                                     tags$hr(),
                                                                                     fluidRow(column(3, downloadButton("down_peplfq", "Download results")),
                                                                                              column(3, selectInput("format_peplfq", "Select a format",
                                                                                                                    choices = c("txt", "csv", "xlsx"),
                                                                                                                    selected = "txt"))
                                                                                              )
                                                                                     )
                                                                    )
                                                                )
                                                       )
                                      ),
                             tabPanel("Protein group and genes",
                                      conditionalPanel(condition = "!output.reportdata_up",
                                                       h2("Import a report file in the tab 'Import your data' first !")
                                                       ),
                                      conditionalPanel(condition = "output.reportdata_up",
                                                       tags$hr(),
                                                       fluidRow(box(title = "Protein group", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                                    tags$u(h3("Get your protein group file (will use the MaxLFQ algorithm)")),
                                                                    tags$hr(),

                                                                    fluidRow(column(3, numericInput("qv_pg", "Choose the q-value to filter the proteins",
                                                                                                    min = 0, max = 1, step = 0.01, value = 0.01)),
                                                                             column(3, numericInput("qvprot_pg", "Choose the protein.q-value to filter the proteins",
                                                                                                    min = 0, max = 1, step = 0.01, value = 1)),
                                                                             column(3, numericInput("qvpg_pg", "Choose the protein-group q-value to filter the proteins",
                                                                                                    min = 0, max = 1, step = 0.01, value = 0.01)),
                                                                             column(3, numericInput("qvgg_pg", "Choose the gene-groupe q-value to filter the proteins",
                                                                                                    min = 0, max = 1, step = 0.01, value = 1))
                                                                             ),
                                                                    radioButtons("wLFQ_pg", "",
                                                                                 choices = c("Use fast MaxLFQ from iq package (log2 transformed)" = "iq",
                                                                                             "Use MaxLFQ from diann package" = "diann"),
                                                                                 selected = "iq",
                                                                                 inline = TRUE),
                                                                    fluidRow(column(3, checkboxInput("onlycountall_pg", "Only keep peptides counts all", TRUE)),
                                                                             column(3, checkboxInput("protypiconly_pg", "Proteotypic only", TRUE)),
                                                                             column(3, checkboxInput("Top3_pg", "Get Top3 quantification", TRUE)),
                                                                             column(3, checkboxInput("iBAQ_pg", "Get iBAQ quantification", TRUE))
                                                                             ),
                                                                    conditionalPanel(condition = "input.iBAQ_pg",
                                                                                     fluidRow(column(4, checkboxInput("fasta_pg", "Import you own FASTA files; if not, search on swissprot.", TRUE),
                                                                                                     conditionalPanel(condition = "input.fasta_pg",
                                                                                                                      fileInput("fastafile_pg", "Import your FATSA files", multiple = TRUE),
                                                                                                                      ),
                                                                                                     conditionalPanel(condition = "!input.fasta_pg",
                                                                                                                      selectizeInput("species_pg", "Choose a species",
                                                                                                                                     choices = NULL)
                                                                                                                      )
                                                                                                     ),
                                                                                              column(4, sliderInput("peplen_pg", "Choose the min and max peptide length",
                                                                                                                    min = 0, max = 100, value = c(5,36), step = 1)
                                                                                                     ),
                                                                                              column(4, selectInput("enzyme_pg", "Choose an enzyme", choices = c("trypsin", "lys-c"), selected = "trypsin")
                                                                                                     )
                                                                                              ),
                                                                                     tags$hr(),
                                                                                     textOutput("diag_getseq"),
                                                                                     tags$hr(),
                                                                                     ),
                                                                    fluidRow(column(3, selectInput("modif_pg", "Choose some modifications to remove (if none selected, no filtering)",
                                                                                                   choices = "",
                                                                                                   multiple = TRUE))
                                                                             ),
                                                                    actionButton("go_pg", "Start calculation", class = "btn-primary"),
                                                                    tags$hr(),
                                                                    textOutput("info_pg"),
                                                                    conditionalPanel(condition = "output.proteins_up",
                                                                                     DT::dataTableOutput("res_pg"),
                                                                                     tags$hr(),
                                                                                     fluidRow(column(3, downloadButton("down_pg", "Download results")),
                                                                                              column(3, selectInput("format_pg", "Select a format",
                                                                                                                    choices = c("txt", "csv", "xlsx"),
                                                                                                                    selected = "txt"))
                                                                                              )
                                                                                     )
                                                                    ),
                                                                shinyjs::useShinyjs()
                                                                ),
                                                       fluidRow(box(title = "Unique genes", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                                    tags$u(h3("Get your unique genes file")),
                                                                    tags$hr(),

                                                                    fluidRow(column(3, numericInput("qv_gg", "Choose the q-value to filter the genes",
                                                                                                    min = 0, max = 1, step = 0.01, value = 0.01)),
                                                                             column(3, numericInput("qvprot_gg", "Choose the protein.q-value to filter the genes",
                                                                                                    min = 0, max = 1, step = 0.01, value = 1)),
                                                                             column(3, numericInput("qvpg_gg", "Choose the protein-group q-value to filter the genes",
                                                                                                    min = 0, max = 1, step = 0.01, value = 0.01)),
                                                                             column(3, numericInput("qvgg_gg", "Choose the gene-groupe q-value to filter the genes",
                                                                                                    min = 0, max = 1, step = 0.01, value = 1))
                                                                             ),
                                                                    fluidRow(column(3, checkboxInput("onlycountall_gg", "Only keep peptides counts all", TRUE)),
                                                                             column(3, checkboxInput("protypiconly_gg", "Proteotypic only", TRUE)),
                                                                             column(3, checkboxInput("Top3_gg", "Get Top3 quantification", TRUE))
                                                                             ),
                                                                    fluidRow(column(3, selectInput("modif_gg", "Choose some modifications to remove (if none selected, no filtering)",
                                                                                                   choices = "",
                                                                                                   multiple = TRUE))
                                                                             ),
                                                                    actionButton("go_gg", "Start calculation", class = "btn-primary"),
                                                                    tags$hr(),
                                                                    textOutput("info_gg"),
                                                                    conditionalPanel(condition = "output.genes_up",
                                                                                     DT::dataTableOutput("res_gg"),
                                                                                     tags$hr(),
                                                                                     fluidRow(column(3, downloadButton("down_gg", "Download results")),
                                                                                              column(3, selectInput("format_gg", "Select a format",
                                                                                                                    choices = c("txt", "csv", "xlsx"),
                                                                                                                    selected = "txt"))
                                                                                              )
                                                                                     )
                                                                    )
                                                                )
                                                       )
                                      ),
                             tabPanel("Data visualization",
                                      tags$hr(),
                                      tabsetPanel(type = "tabs",
                                                  tabPanel("Check results",
                                                           tags$hr(),
                                                            fluidRow(box(title = "Data choice", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                                         conditionalPanel(condition = "output.reportdata_up",
                                                                                          radioButtons("choice_visu", "",
                                                                                                       choices = c("Visualize imported data" = "base",
                                                                                                                   "Import your own data" = "dat"),
                                                                                                       selected = "dat",
                                                                                                       inline = TRUE),
                                                                                          conditionalPanel(condition = "input.choice_visu == 'base'",
                                                                                                           radioButtons("bdata_visu", "",
                                                                                                                        choices = c("Protein group",
                                                                                                                                    "Unique genes",
                                                                                                                                    "Peptides",
                                                                                                                                    "Peptides.MaxLFQ",
                                                                                                                                    "Precursors"),
                                                                                                                        selected = "Protein group",
                                                                                                                        inline = TRUE)
                                                                                                           )
                                                                                          ),
                                                                         conditionalPanel(condition = "!output.reportdata_up | input.choice_visu == 'dat'",
                                                                                          fileInput("ydata_visu", "Import your data. The protein group file, for example.")
                                                                                          )
                                                                         )
                                                                     ),

                                                            conditionalPanel(condition = "output.visudata_up",
                                                                             fluidRow(box(title = "Visualiation", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                                                          conditionalPanel(condition = "output.top3_or_ibaq",
                                                                                                           selectInput("dtype_visu", "Choose the type of data you want to visualize.",
                                                                                                                       choices = c("intensity"))
                                                                                                           ),
                                                                                          tabsetPanel(type = "tabs",
                                                                                                      tabPanel("Heatmap",
                                                                                                               tags$hr(),
                                                                                                               fluidRow(column(3, selectInput("transfo_visu", "Choose a data transformation",
                                                                                                                                              choices = c("Log2" = "log2",
                                                                                                                                                          "Z-score on the proteins" = "z.score_protein",
                                                                                                                                                          "Z-score on the fractions" = "z.score_fraction",
                                                                                                                                                          "None" = "none"), selected = "none")),
                                                                                                                        column(3, checkboxInput("prval_visu", "Print values on blocks", FALSE),
                                                                                                                               sliderInput("maxna", "Choose the maximum number of missing values per rows",
                                                                                                                                           value = 0, min = 0, step = 1, max = 3)),
                                                                                                                        conditionalPanel(condition = "!output.reportdata_up | input.choice_visu == 'dat'",
                                                                                                                                         column(3, textInput("nmid_visu", "Type the name of the column that contains
                                                                                                                                        the IDs"))
                                                                                                                                         ),
                                                                                                                        column(3, actionButton("seeheat_visu", "See heatmap", class = "btn-primary"))
                                                                                                                        ),
                                                                                                               tags$hr(),
                                                                                                               withSpinner(plotlyOutput("heat_visu", height = "800px"), type = 6)
                                                                                                               ),
                                                                                                      tabPanel("Density",
                                                                                                               tags$hr(),
                                                                                                               fluidRow(column(3, selectInput("transfoD_visu", "Choose a data transformation",
                                                                                                                                              choices = c("Log2" = "log2",
                                                                                                                                                          "None" = "none"), selected = "none")),
                                                                                                                        column(3, checkboxInput("area_visu", "Print area of the density curves", TRUE)),
                                                                                                                        column(3, textInput("titD_visu", "Choose a title for your plot (can be NULL)")),
                                                                                                                        column(3, actionButton("seedens_visu", "See density plot", class = "btn-primary"))
                                                                                                               ),
                                                                                                               tags$hr(),
                                                                                                               withSpinner(plotOutput("dens_visu", height = "800px"), type = 6),
                                                                                                               downloadButton("down_dens", "Download density plot")
                                                                                                               ),
                                                                                                      tabPanel("MDS",
                                                                                                               tags$hr(),
                                                                                                               fluidRow(column(4, selectInput("transfoM_visu", "Choose a data transformation",
                                                                                                                                              choices = c("Log2" = "log2",
                                                                                                                                                          "None" = "none"), selected = "none")),
                                                                                                                        column(4, textInput("titM_visu", "Choose a title for your plot (can be NULL)")),
                                                                                                                        column(4, actionButton("seemds_visu", "See MDS plot", class = "btn-primary"))
                                                                                                               ),
                                                                                                               tags$hr(),
                                                                                                               withSpinner(plotOutput("mds_visu", height = "800px"), type = 6),
                                                                                                               downloadButton("down_mds", "Download MDS plot")
                                                                                                               ),
                                                                                                      tabPanel("Proportion of valid values",
                                                                                                               tags$hr(),
                                                                                                               fluidRow(column(3, selectInput("transfoPVV_visu", "Choose a data transformation",
                                                                                                                                              choices = c("Log2" = "log2",
                                                                                                                                                          "None" = "none"), selected = "none"),
                                                                                                                               numericInput("propcutPVV_visu", "Choose the minimum proportion to show on the graph",
                                                                                                                                           value = 0, min = 0, max = 1, step = 0.05)),
                                                                                                                        column(3, textInput("titPVV_visu", "Choose a title for your plot (can be NULL)")),
                                                                                                                        column(3, textInput("designPVV_visu", "Type the design of your experiment that match your columns names,
                                                                                                                                                               separated by a comma.
                                                                                                                                                               For example, if you have columns named '37_B1_Treatment', type 'temperature,replicate,condition'.")),
                                                                                                                        column(3, selectInput("checkPVV_visu", "Choose the grouping variable (same as in your design; like 'replicate' for example)",
                                                                                                                                              choices = NULL))
                                                                                                                        ),
                                                                                                               fluidRow(column(3, actionButton("seePVV_visu", "See PVV plot", class = "btn-primary"))),
                                                                                                               tags$hr(),
                                                                                                               withSpinner(plotOutput("PVV_visu", height = "800px"), type = 6),
                                                                                                               downloadButton("down_PVV", "Download Proportion valid values plot")
                                                                                                               ),
                                                                                                      tabPanel("Retention time",
                                                                                                               conditionalPanel(condition = "input.choice_visu == 'base'",
                                                                                                                                tags$hr(),
                                                                                                                                fluidRow(column(4, actionButton("seert_visu", "See RT vs iRT", class = "btn-primary"))
                                                                                                                                         ),
                                                                                                                                tags$hr(),
                                                                                                                                withSpinner(plotlyOutput("rt1_visu", height = "800px"), type = 6),
                                                                                                                                withSpinner(plotlyOutput("rt2_visu", height = "800px"), type = 6)
                                                                                                                                ),
                                                                                                               conditionalPanel(condition = "input.choice_visu == 'dat'",
                                                                                                                                h3("You need to import the report file from DIA nn to use this tab.
                                                                                                                                    For this, go back to the first tab 'Import your data'"))
                                                                                                               ),
                                                                                                      tabPanel("Proteotypic",
                                                                                                               conditionalPanel(condition = "input.choice_visu == 'base'",
                                                                                                                                tags$hr(),
                                                                                                                                fluidRow(column(4, actionButton("seeptyp_visu", "See proteotypic proportion", class = "btn-primary"))
                                                                                                                                ),
                                                                                                                                tags$hr(),
                                                                                                                                withSpinner(plotOutput("ptyp1_visu", height = "600px"), type = 6),
                                                                                                                                withSpinner(plotOutput("ptyp2_visu", height = "600px"), type = 6)
                                                                                                               ),
                                                                                                               conditionalPanel(condition = "input.choice_visu == 'dat'",
                                                                                                                                h3("You need to import the report file from DIA nn to use this tab.
                                                                                                                                    For this, go back to the first tab 'Import your data'"))
                                                                                                               )
                                                                                                      )
                                                                                          )
                                                                                      )
                                                                             )
                                                            ),
                                                  tabPanel("Check other reports",
                                                           tags$hr(),
                                                           fluidRow(box(title = "Window selection", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
                                                                        fluidRow(column(6, shinyFilesButton("replib_wsel", label = "Select your report file", title = "Please select a file",
                                                                                                            icon = icon("file"),
                                                                                                            multiple = TRUE, viewtype = "detail", buttonType = "primary", class = "btn-lg"))
                                                                        ),
                                                                        conditionalPanel(condition = "output.libdata_up",
                                                                                         tags$hr(),
                                                                                         fluidRow(column(4, numericInput("bins_wsel", "Select the number of windows you want to have",
                                                                                                                         25, min = 5, step = 1)),
                                                                                                  column(4, checkboxInput("perfrac_wsel", "Average the best windows from each fraction", FALSE))
                                                                                         ),
                                                                                         fluidRow(column(4, actionButton("go_wsel", "Get best windows", class = "btn-primary")),
                                                                                                  column(6, radioButtons("whichplot_wsel", "", choices = c("Visualize global distribution" = "all",
                                                                                                                                                           "Visualize distribution from each fraction" = "frac"),
                                                                                                                         selected = "all",
                                                                                                                         inline = TRUE)),
                                                                                                  ),
                                                                                         tags$hr(),
                                                                                         textOutput("bstw_wsel"),
                                                                                         tags$hr(),
                                                                                         conditionalPanel(condition = "input.whichplot_wsel == 'all'",
                                                                                                          withSpinner(plotOutput("hist_wsel", height = "800px"), type = 6),
                                                                                                          downloadButton("downhist_wsel", "Download plot")),
                                                                                         conditionalPanel(condition = "input.whichplot_wsel == 'frac'",
                                                                                                          withSpinner(plotOutput("hist_wsel_frac", height = "800px"), type = 6),
                                                                                                          downloadButton("downhist_wsel_frac", "Download plot")),
                                                                                         shinyjs::useShinyjs()
                                                                                         )
                                                                        )
                                                                    )
                                                           )
                                                  )
                                      )
                             )
                 )
               )
             ),

    bslib::nav_item(a(href = "mailto:marco.gerault@gmail.com",
                      icon("envelope"),
                      title = "Any questions, suggestions or bug report ? Feel free to send me an e-mail !")
                    ),
    bslib::nav_item(a(href = "https://youtu.be/vfvh15Q93eU",
                      icon("question-circle"),
                      title = "See the tutorial video")
                    )
    )
)

server <- function(input, output, session){
  observe({
    updateSelectizeInput(session, "species_pg", choices = DIAgui::all_species, selected = "HOMO SAPIENS", server = TRUE)
  })

  ### REPORT FILE
  pth <- str_split(WD, "/")[[1]]
  pth <- paste(pth[1:4][!is.na(pth[1:4])], collapse = "/")
  names(pth) <- pth
  volumes <- c(Home = WD, "R Installation" = R.home(), getVolumes()(), pth)
  shinyFileChoose(input, "rep_tsv", roots = volumes, session = session)

  output$filepaths <- renderPrint({
    if (is.integer(input$rep_tsv)) {
      cat("No files have been selected (shinyFileChoose)")
    } else {
      parseFilePaths(volumes, input$rep_tsv)
    }
  })

  report_data <- reactive({
    if (is.integer(input$rep_tsv)) {
      return(NULL)
    }
    else {
      File <- parseFilePaths(volumes, input$rep_tsv)
    }
    File <- parseFilePaths(volumes, input$rep_tsv)
    if(is.null(File))
      return(NULL)

    showNotification("Getting your data, this may take a while.", type = "message")
    diann_load(File$datapath)
  })
  Report_data <- reactiveValues(
    d = NULL
  )
  observe({
    Report_data$d <- report_data()
  })
  output$reportdata_up <- reactive({
    return(!is.null(Report_data$d))
  })
  outputOptions(output, "reportdata_up", suspendWhenHidden = FALSE)

  output$df_report <- DT::renderDataTable({
    DT::datatable(Report_data$d,
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    htmltools::strong('Report from DIA-nn')
                    ),
                  rownames = FALSE,
                  options = list(lenghtMenu = c(10,20,30), pageLength = 10,
                                 scrollX = TRUE)
                  )
    })


  output$frac_dat <-  renderUI({
    if(!is.null(Report_data$d)){
      fr <- unique(Report_data$d$File.Name)
    }
    else{
      fr <- unique(report_data()$File.Name)
    }
    if(length(fr) == 1){
      HTML(paste0("<h3><u>The current fraction name is :</u> ", fr, ".</h3>"
                  )
           )
    }
    else{
      HTML(paste0("<h3><u>The current fraction names are :</u> ", paste(paste(fr[1:(length(fr)-1)], collapse = ", "),
                                                                        "and", fr[length(fr)]), ".</h3>")
           )
    }

  })

  observeEvent(input$change_dat, {
    if(!is.null(Report_data$d)){
      report <- Report_data$d
      cd <- unique(report$File.Name)
    }

    if(input$chorename_dat == "New names"){
      if(str_length(str_remove_all(input$newfrac_dat, " ")) == 0){
        showNotification("Type something !", type = "error", duration = 2)
      }
      else if(str_count(input$newfrac_dat, ",") != (length(cd) - 1)){
        showNotification("You need to type the same number of new names as you have fractions,
                        even if it's empty", type = "error", duration = 4)
      }
      else{
        showNotification("Checking new names", type = "message", duration = 2)
        nm <- input$newfrac_dat
        nm <- str_split(nm, ",")[[1]]
        nm <- str_remove_all(nm, " ")
        if(sum(str_detect(nm, "/")) > 0){
          showNotification("The character '/' is not alllowed. Please, verify your new names.", type = "error")
        }
        else{
          showNotification("Start changing names", type = "message", duration = 3)
          for(i in 1:length(cd)){
            if(str_length(nm[i]) == 0){
              nm[i] <- cd[i]
            }
          }
          change <- cd[!(nm %in% cd)]
          if(!purrr::is_empty(change) & !is.null(change)){
            new <- nm[!(nm %in% cd)]
            showNotification(paste("You decided to change :", paste(change, collapse = ", "),
                                   "In :", paste(new, collapse = ", ")), type = "message", duration = 7)

            for(i in 1:length(change)){
              report$File.Name[which(report$File.Name == change[i])] <- new[i]
            }
            Report_data$d <- report

            showNotification("Names changed !", type = "message", duration = 3)
          }
          else{
            showNotification("You didn't make any changement !", type = "error")
          }
        }
      }
    }
    else if(input$chorename_dat == "Remove path"){
      showNotification("Start changing names", type = "message", duration = 3)
      if(input$whattorm_dat == 1){
        nm <- lapply(cd, function(x) as.data.frame(t(str_split_fixed(x, "", str_length(x)))))
        N <- max(as.numeric(lapply(nm, nrow)))
        nm <- lapply(nm, function(x) {x <- rbind(x, t(t(rep("", N - nrow(x))))); x})
        nm <- do.call(cbind, nm)
        names(nm) <- paste0("F", 1:length(nm))
        nm$keep <- apply(nm, 1, function(x) length(unique(x)) != 1)
        nm <- as.list(nm[nm$keep, -ncol(nm)])
        nm <- lapply(nm, function(x) paste(x, collapse = ""))

        for(i in 1:length(nm)){
          report$File.Name[which(report$File.Name == cd[i])] <- nm[[i]]
        }
      }
      else if(input$whattorm_dat == 2){
        report$File.Name <- str_remove_all(report$File.Name, ".{1,}\\\\")
      }
      else if(input$whattorm_dat == 3){
        report$File.Name <- str_remove_all(report$File.Name, ".{1,}\\\\")
        report$File.Name <- str_remove_all(report$File.Name, "\\..{1,}")
      }
      Report_data$d <- report

      showNotification("Names changed !", type = "message", duration = 3)
    }

  })
    ### PRECURSORS
    observe({
      if(!is.null(Report_data$d)){
        modif <- Report_data$d$Modified.Sequence
        modif <- unique(stringr::str_extract(modif, "(?<=\\().+?(?=\\))"))
        modif <- modif[!is.na(modif)]

        updateSelectInput(session, "modif_prec", choices = modif)
      }
    })
    precu_ev <- reactiveValues(
      x = NULL
    )
    precu <- reactive({
      df <- Report_data$d
      if(!is.null(input$modif_prec)){
        idx_modif <- stringr::str_which(df$Modified.Sequence, paste(input$modif_prec, collapse = "|"))
        if(input$remmodif_prec == "rem"){
          df <- df[-idx_modif,]
        }
        else if(input$remmodif_prec == "keep"){
          df <- df[idx_modif,]
        }
      }
      d <- diann_matrix(df,
                        proteotypic.only = input$protypiconly_prec,
                        q = input$qv_prec,
                        protein.q = input$qvprot_prec,
                        pg.q = input$qvpg_prec,
                        gg.q = input$qvgg_prec,
                        method = "max")
    })

    observeEvent(input$go_prec, {
      withCallingHandlers({
        shinyjs::html("info_prec", "")
        message("Calculation...")
        showNotification("Getting the precursors tab", type = "message", duration = 2)
        precu_ev$x <- precu()
        message("Done !")
      },
      message = function(m) {
        shinyjs::html(id = "info_prec", html = paste(m$message, "<br>", sep = ""), add = FALSE)
      }
      )
    })
    output$precursor_up <- reactive({
      return(!is.null(precu_ev$x))
    })
    outputOptions(output, "precursor_up", suspendWhenHidden = FALSE)


    output$res_prec <- DT::renderDataTable({
      DT::datatable(precu_ev$x,
                    caption = htmltools::tags$caption(
                      style = 'caption-side: top; text-align: left;',
                      htmltools::strong('Precursors')
                    ),
                    rownames = FALSE,
                    options = list(lenghtMenu = c(10,20,30), pageLength = 10,
                                   scrollX = TRUE)
      )
    })
    output$down_prec <- downloadHandler(
      filename = function() {
        paste0("Precursors_dia_", Sys.Date(), ".", input$format_prec)
      },
      content = function(file){
        if(input$format_prec == "xlsx"){
          openxlsx::write.xlsx(precu_ev$x, file, rowNames = FALSE)
        }
        else if(input$format_prec == "csv"){
          write.csv(precu_ev$x, file, row.names =  FALSE, quote = FALSE)
        }
        else if(input$format_prec == "txt"){
          write.table(precu_ev$x, file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
      }
    )

    ### PEPTIDE
    observe({
      if(!is.null(Report_data$d)){
        modif <- Report_data$d$Modified.Sequence
        modif <- unique(stringr::str_extract(modif, "(?<=\\().+?(?=\\))"))
        modif <- modif[!is.na(modif)]

        updateSelectInput(session, "modif_pep", choices = modif)
      }
    })
    pep_ev <- reactiveValues(
      x = NULL
    )
    pep <- reactive({
      df <- Report_data$d
      if(!is.null(input$modif_pep)){
        idx_modif <- stringr::str_which(df$Modified.Sequence, paste(input$modif_pep, collapse = "|"))
        if(input$remmodif_pep == "rem"){
          df <- df[-idx_modif,]
        }
        else if(input$remmodif_pep == "keep"){
          df <- df[idx_modif,]
        }
      }
      d <- diann_matrix(df,
                        id.header = input$centercol_pep,
                        proteotypic.only = input$protypiconly_pep,
                        q = input$qv_pep,
                        protein.q = input$qvprot_pep,
                        pg.q = input$qvpg_pep,
                        gg.q = input$qvgg_pep,
                        method = "max")

      if(input$getPTM_pep){
        if(all(c("PTM.Q.Value", "PTM.Site.Confidence","Lib.PTM.Site.Confidence") %in% colnames(df))){
          ptm <- df[,c(colnames(d)[1:5],
                      "PTM.Q.Value", "PTM.Site.Confidence",
                      "Lib.PTM.Site.Confidence")]

          ptm <- ptm %>%
            group_by(Modified.Sequence, Stripped.Sequence,
                     Protein.Group, Protein.Names, Genes) %>%
            summarize(PTM.Q.Value = min(PTM.Q.Value),
                      PTM.Site.Confidence = max(PTM.Site.Confidence),
                      Lib.PTM.Site.Confidence = max(Lib.PTM.Site.Confidence)) %>%
            right_join(d, by = colnames(d)[1:5])

          d <- ptm
        }
        else{
          showNotification("Your report doesn't contain the PTM informations !", type = "warning", duration = 5)
        }
      }

      d
    })
    observeEvent(input$go_pep, {
      withCallingHandlers({
        shinyjs::html("info_pep", "")
        if(input$getPTM_pep & input$centercol_pep == "Stripped.Sequence"){
          showNotification("If you want to extract the PTM q.values,
                            you should select the 'Modified.Sequence'", type = "error", duration = 5)
        }
        else{
          message("Calculation...")
          showNotification("Getting the peptides tab", type = "message", duration = 2)
          pep_ev$x <- pep()
          message("Done !")
        }
      },
      message = function(m) {
        shinyjs::html(id = "info_pep", html = paste(m$message, "<br>", sep = ""), add = FALSE)
      }
      )

    })
    output$peptide_up <- reactive({
      return(!is.null(pep_ev$x))
    })
    outputOptions(output, "peptide_up", suspendWhenHidden = FALSE)

    output$res_pep <- DT::renderDataTable({
      DT::datatable(pep_ev$x,
                    caption = htmltools::tags$caption(
                      style = 'caption-side: top; text-align: left;',
                      htmltools::strong('Peptides')
                    ),
                    rownames = FALSE,
                    options = list(lenghtMenu = c(10,20,30), pageLength = 10,
                                   scrollX = TRUE)
      )
    })
    output$down_pep <- downloadHandler(
      filename = function() {
        paste0("Peptides_dia_", Sys.Date(), ".", input$format_pep)
      },
      content = function(file){
        if(input$format_pep == "xlsx"){
          openxlsx::write.xlsx(pep_ev$x, file, rowNames = FALSE)
        }
        else if(input$format_pep == "csv"){
          write.csv(pep_ev$x, file, row.names =  FALSE, quote = FALSE)
        }
        else if(input$format_pep == "txt"){
          write.table(pep_ev$x, file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
      }
    )

    ### PEPTIDE MaxLFQ
    observe({
      if(!is.null(Report_data$d)){
        modif <- Report_data$d$Modified.Sequence
        modif <- unique(stringr::str_extract(modif, "(?<=\\().+?(?=\\))"))
        modif <- modif[!is.na(modif)]

        updateSelectInput(session, "modif_peplfq", choices = modif)
      }
    })
    peplfq_ev <- reactiveValues(
      x = NULL
    )
    peplfq <- reactive({
      df <- Report_data$d
      if(!is.null(input$modif_peplfq)){
        idx_modif <- stringr::str_which(df$Modified.Sequence, paste(input$modif_peplfq, collapse = "|"))
        if(input$remmodif_peplfq == "rem"){
          df <- df[-idx_modif,]
        }
        else if(input$remmodif_peplfq == "keep"){
          df <- df[idx_modif,]
        }
      }
      if(input$protypiconly_peplfq){
        df <- df[which(df[["Proteotypic"]] != 0), ]
      }
      df <- df %>% dplyr::filter(Q.Value <= input$qv_peplfq & PG.Q.Value <= input$qvpg_peplfq & Protein.Q.Value <= input$qvprot_peplfq & GG.Q.Value <= input$qvgg_peplfq)
      if(input$wLFQ_peplfq == "diann"){
        d <- diann_maxlfq(df,
                          group.header = input$centercol_peplfq,
                          id.header = "Precursor.Id",
                          quantity.header = "Precursor.Normalised",
                          count_pep = FALSE
        )
      }
      else if(input$wLFQ_peplfq == "iq"){
        d <- iq::preprocess(df,
                            intensity_col = "Precursor.Normalised",
                            primary_id = input$centercol_peplfq,
                            sample_id  = "File.Name",
                            secondary_id = "Precursor.Id",
                            median_normalization = FALSE,
                            pdf_out = NULL)
        d <- iq::fast_MaxLFQ(d)
        d <- d$estimate
        d <- as.data.frame(d)
      }
      nc <- ncol(d)
      d <- d[order(rownames(d)), , drop=FALSE]

      df <- df[(df$Modified.Sequence %in% rownames(d)),]
      df <- df[order(df$Modified.Sequence),]
      m <- unique(df[,c("Modified.Sequence", "Stripped.Sequence", "Protein.Group", "Protein.Names", "Genes")])
      m <- m[order(m$Stripped.Sequence),]

      d <- cbind(d,m)
      rownames(d) <- 1:nrow(d)
      d <- d[,c((nc+1):ncol(d), 1:nc)]

      if(input$getPTM_peplfq){
        if(all(c("PTM.Q.Value", "PTM.Site.Confidence","Lib.PTM.Site.Confidence") %in% colnames(df))){
          ptm <- df[,c(colnames(d)[1:5],
                       "PTM.Q.Value", "PTM.Site.Confidence",
                       "Lib.PTM.Site.Confidence")]

          ptm <- ptm %>%
            group_by(Modified.Sequence, Stripped.Sequence,
                     Protein.Group, Protein.Names, Genes) %>%
            summarize(PTM.Q.Value = min(PTM.Q.Value),
                      PTM.Site.Confidence = max(PTM.Site.Confidence),
                      Lib.PTM.Site.Confidence = max(Lib.PTM.Site.Confidence)) %>%
            right_join(d, by = colnames(d)[1:5])

          d <- ptm
        }
        else{
          showNotification("Your report doesn't contain the PTM informations !", type = "warning", duration = 5)
        }
      }

      d
    })
    observeEvent(input$go_peplfq, {
      if(input$centercol_peplfq == "Stripped.Sequence"){
        showModal(
          modalDialog(
            title = "Using 'Stripped.Sequence' is not advised when using MaxLFQ;
                     number of rows may not match so it might crash.
                     If so, use 'Modified.Sequence' instead.",
            footer = tagList(actionButton("confirm_peplfq", "Confirm"),
                             modalButton("Cancel"))
          )
        )
      }
      else{
        withCallingHandlers({
          shinyjs::html("info_peplfq", "")
          if(input$getPTM_peplfq & input$centercol_peplfq == "Stripped.Sequence"){
            showNotification("If you want to extract the PTM q.values,
                            you should select the 'Modified.Sequence'", type = "error", duration = 5)
          }
          else{
            message("Calculation...")
            showNotification(paste("Getting the peptides tab using the MaxLFQ algorithm from", input$wLFQ_peplfq, "package"), type = "message")
            peplfq_ev$x <- peplfq()
            message("Done !")
          }
        },
        message = function(m) {
          shinyjs::html(id = "info_peplfq", html = paste(m$message, "<br>", sep = ""), add = FALSE)
        }
        )
      }
    })
    observeEvent(input$confirm_peplfq, {
      removeModal()

      withCallingHandlers({
        shinyjs::html("info_peplfq", "")
        message("Calculation...")
        showNotification(paste("Getting the peptides tab using the MaxLFQ algorithm from", input$wLFQ_peplfq, "package"), type = "message")
        peplfq_ev$x <- peplfq()
        message("Done !")
      },
      message = function(m) {
        shinyjs::html(id = "info_peplfq", html = paste(m$message, "<br>", sep = ""), add = FALSE)
      }
      )

    })
    output$peptideLFQ_up <- reactive({
      return(!is.null(peplfq_ev$x))
    })
    outputOptions(output, "peptideLFQ_up", suspendWhenHidden = FALSE)

    output$res_peplfq <- DT::renderDataTable({
      DT::datatable(peplfq_ev$x,
                    caption = htmltools::tags$caption(
                      style = 'caption-side: top; text-align: left;',
                      htmltools::strong('Peptides, using the MaxLFQ algorithm')
                    ),
                    rownames = FALSE,
                    options = list(lenghtMenu = c(10,20,30), pageLength = 10,
                                   scrollX = TRUE)
      )
    })
    output$down_peplfq <- downloadHandler(
      filename = function() {
        paste0("PeptidesMaxLFQ_dia_", Sys.Date(), ".", input$format_peplfq)
      },
      content = function(file){
        if(input$format_peplfq == "xlsx"){
          openxlsx::write.xlsx(peplfq_ev$x, file, rowNames = FALSE)
        }
        else if(input$format_peplfq == "csv"){
          write.csv(peplfq_ev$x, file, row.names =  FALSE, quote = FALSE)
        }
        else if(input$format_peplfq == "txt"){
          write.table(peplfq_ev$x, file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
      }
    )

    ### PROTEINS
    observe({
      if(!is.null(Report_data$d)){
        modif <- Report_data$d$Modified.Sequence
        modif <- unique(stringr::str_extract(modif, "(?<=\\().+?(?=\\))"))
        modif <- modif[!is.na(modif)]

        updateSelectInput(session, "modif_pg", choices = modif)
      }
    })
    fasta <- visu_data <- reactive({
      fts <- NULL
      if(input$fasta_pg){
        File <- input$fastafile_pg
        if(is.null(File))
          return(NULL)

        fts <- File$datapath
      }
      else{
        fts <- NULL
      }
      fts
    })


    pg_ev <- reactiveValues(
      x = NULL
    )
    pg <- reactive({
      withCallingHandlers({
        shinyjs::html("diag_getseq", "")

      df <- Report_data$d
      if(!is.null(input$modif_pg)){
        idx_modif <- stringr::str_which(df$Modified.Sequence, paste(input$modif_pg, collapse = "|"))
        df <- df[-idx_modif,]
      }
      brut <- df
      if(input$protypiconly_pg){
        df <- df[which(df[["Proteotypic"]] != 0), ]
      }
      df <- df %>% dplyr::filter(Q.Value <= input$qv_pg & PG.Q.Value <= input$qvpg_pg & Protein.Q.Value <= input$qvprot_pg & GG.Q.Value <= input$qvgg_pg)
      n_cond <- length(unique(df$File.Name))

      if(input$wLFQ_pg == "diann"){
        d <- diann_maxlfq(df,
                          group.header="Protein.Group",
                          id.header = "Precursor.Id",
                          quantity.header = "Precursor.Normalised",
                          only_countsall = input$onlycountall_pg,
                          Top3 = input$Top3_pg
                          )
      }
      else if(input$wLFQ_pg == "iq"){
        d <- iq::preprocess(df,
                            intensity_col = "Precursor.Normalised",
                            primary_id = "Protein.Group",
                            sample_id  = "File.Name",
                            secondary_id = "Precursor.Id",
                            median_normalization = FALSE,
                            pdf_out = NULL)

        pc <- d %>% dplyr::group_by(protein_list, sample_list) %>%
          dplyr::mutate("countpep" = length(unique(id)))
        pc <- unique(pc[,c("protein_list", "sample_list", "countpep")])
        pc <- tidyr::spread(pc, sample_list, countpep)
        pc[is.na(pc)] <- 0
        pc <- as.data.frame(pc)
        rownames(pc) <- pc$protein_list
        pc$protein_list <- NULL
        pc <- pc[order(rownames(pc)), , drop=FALSE]
        colnames(pc) <- paste0("pep_count_", colnames(pc))
        pc$peptides_counts_all <- unname(apply(pc, 1, max))
        pc <- pc[,c(ncol(pc), 1:(ncol(pc)-1))]

        if(input$Top3_pg){
          Top3 <- iq::preprocess(df,
                              intensity_col = "Precursor.Quantity",
                              primary_id = "Protein.Group",
                              sample_id  = "File.Name",
                              secondary_id = "Precursor.Id",
                              median_normalization = FALSE,
                              pdf_out = NULL)
          Top3 <- Top3 %>% dplyr::group_by(protein_list) %>%
            tidyr::spread(sample_list, quant)
          Top3$id <- NULL
          top3_f <- function(x){
            if(sum(!is.na(x)) < 3){
              x <- NA
            }
            else{
              x <- x[order(x, decreasing = TRUE)][1:3]
              x <- mean(x)
            }
          }
          Top3 <- Top3 %>% dplyr::group_by(protein_list) %>%
            dplyr::summarise(dplyr::across(dplyr::everything(), top3_f))

          Top3 <- as.data.frame(Top3)
          rownames(Top3) <- Top3$protein_list
          Top3$protein_list <- NULL
          colnames(Top3) <- paste0("Top3_", colnames(Top3))
          Top3 <- Top3[order(rownames(Top3)), , drop=FALSE]
        }

        d <- iq::fast_MaxLFQ(d)
        d <- d$estimate
        d <- as.data.frame(d)
        d <- d[order(rownames(d)), , drop=FALSE]
        if(input$Top3_pg){
          d <- cbind(d, Top3)
        }
        if(input$onlycountall_pg){
          d$peptides_counts_all <- pc$peptides_counts_all
        }
        else{
          d <- cbind(d, pc)
        }
      }
      nc <- ncol(d)
      d$Protein.Group <- rownames(d)
      rownames(d) <- 1:nrow(d)
      df <- df[(df$Protein.Group %in% d$Protein.Group),]
      df <- df[order(df$Protein.Group),]
      d$Protein.Names <- unique(df[,c("Protein.Group", "Protein.Names")])$Protein.Names
      d$First.Protein.Description <- unique(df[,c("Protein.Group", "First.Protein.Description")])$First.Protein.Description
      d$Genes <- unique(df[,c("Protein.Group", "Genes")])$Genes
      d <- d[,c((nc+1):ncol(d), 1:nc)]

      if(input$iBAQ_pg){
        if(input$fasta_pg){
          d_seq <- getallseq(pr_id = d$Protein.Group,
                             fasta_file = TRUE,
                             bank_name = fasta())
        }
        else{
          d_seq <- getallseq(pr_id = d$Protein.Group,
                             spec = input$species_pg)
        }

        brut <- diann_matrix(brut, id.header = "Protein.Group",
                             quantity.header = "Precursor.Quantity",
                             proteotypic.only = input$protypiconly_pg,
                             q = input$qv_pg, protein.q = input$qvprot_pg,
                             pg.q = input$qvpg_pg, gg.q = input$qvgg_pg,
                             method = "sum")
        brut$Genes <- NULL
        brut$Protein.Names <- NULL
        brut <- get_iBAQ(brut, proteinDB = d_seq,
                      id_name = "Protein.Group",
                      ecol = 2:(n_cond+1),
                      peptideLength = input$peplen_pg,
                      proteaseRegExp = getProtease(input$enzyme_pg),
                      log2_transformed = input$wLFQ_pg == "iq")
        brut <- brut[,-c(2:(n_cond+1))]
        d <- dplyr::left_join(d, brut, by = "Protein.Group")
      }
      d},
      message = function(m) {
        shinyjs::html(id = "diag_getseq", html = paste(m$message, "<br>", sep = ""), add = FALSE)
      }
      )
    })
    observeEvent(input$go_pg, {
      withCallingHandlers({
        shinyjs::html("info_pg", "")
        message("Calculation...")
        showNotification(paste("Getting the protein group tab using the MaxLFQ algorithm from", input$wLFQ_peplfq, "package"), type = "message")
        pg_ev$x <- pg()
        message("Done !")
      },
      message = function(m) {
        shinyjs::html(id = "info_pg", html = paste(m$message, "<br>", sep = ""), add = FALSE)
      }
      )
    })
    output$proteins_up <- reactive({
      return(!is.null(pg_ev$x))
    })
    outputOptions(output, "proteins_up", suspendWhenHidden = FALSE)

    output$res_pg <- DT::renderDataTable({
      DT::datatable(pg_ev$x,
                    caption = htmltools::tags$caption(
                      style = 'caption-side: top; text-align: left;',
                      htmltools::strong('Protein group')
                      ),
                    rownames = FALSE,
                    options = list(lenghtMenu = c(10,20,30), pageLength = 10,
                                   scrollX = TRUE)
      )
    })
    output$down_pg <- downloadHandler(
      filename = function() {
        paste0("ProteinGroup_dia_", Sys.Date(), ".", input$format_pg)
      },
      content = function(file){
        if(input$format_pg == "xlsx"){
          openxlsx::write.xlsx(pg_ev$x, file, rowNames = FALSE)
        }
        else if(input$format_pg == "csv"){
          write.csv(pg_ev$x, file, row.names =  FALSE, quote = FALSE)
        }
        else if(input$format_pg == "txt"){
          write.table(pg_ev$x, file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
      }
    )

    ### GENES
    observe({
      if(!is.null(Report_data$d)){
        modif <- Report_data$d$Modified.Sequence
        modif <- unique(stringr::str_extract(modif, "(?<=\\().+?(?=\\))"))
        modif <- modif[!is.na(modif)]

        updateSelectInput(session, "modif_gg", choices = modif)
      }
    })
    gg_ev <- reactiveValues(
      x = NULL
    )
    gg <- reactive({
      df <- Report_data$d
      if(!is.null(input$modif_gg)){
        idx_modif <- stringr::str_which(df$Modified.Sequence, paste(input$modif_gg, collapse = "|"))
        df <- df[-idx_modif,]
      }
      d <- diann_matrix(df,
                   id.header="Genes",
                   quantity.header="Genes.MaxLFQ.Unique",
                   proteotypic.only = input$protypiconly_gg,
                   q = input$qv_gg,
                   protein.q = input$qvprot_gg,
                   pg.q = input$qvpg_gg,
                   gg.q = input$qvgg_gg,
                   get_pep = TRUE, only_pepall = input$onlycountall_gg,
                   Top3 = input$Top3_gg,
                   method = "max")
    })
    observeEvent(input$go_gg, {
      withCallingHandlers({
        shinyjs::html("info_gg", "")
        message("Calculation...")
      showNotification("Getting the unique genes tab", type = "message", duration = 2)
      gg_ev$x <- gg()
      message("Done !")
      },
      message = function(m) {
        shinyjs::html(id = "info_gg", html = paste(m$message, "<br>", sep = ""), add = FALSE)
      }
      )
    })
    output$genes_up <- reactive({
      return(!is.null(gg_ev$x))
    })
    outputOptions(output, "genes_up", suspendWhenHidden = FALSE)

    output$res_gg <- DT::renderDataTable({
      DT::datatable(gg_ev$x,
                    caption = htmltools::tags$caption(
                      style = 'caption-side: top; text-align: left;',
                      htmltools::strong('Unique genes')
                    ),
                    rownames = FALSE,
                    options = list(lenghtMenu = c(10,20,30), pageLength = 10,
                                   scrollX = TRUE)
      )
    })
    output$down_gg <- downloadHandler(
      filename = function() {
        paste0("UniqueGenes_dia_", Sys.Date(), ".", input$format_gg)
      },
      content = function(file){
        if(input$format_gg == "xlsx"){
          openxlsx::write.xlsx(gg_ev$x, file, rowNames = FALSE)
        }
        else if(input$format_gg == "csv"){
          write.csv(gg_ev$x, file, row.names =  FALSE, quote = FALSE)
        }
        else if(input$format_gg == "txt"){
          write.table(gg_ev$x, file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
      }
    )


    ### VISUALIZATION
    observe({
      if(!is.null(Report_data$d)){
        updateRadioButtons(session, "choice_visu", selected = "base")
      }
    })
    ## HEATMAP
    visu_data <- reactive({
      df <- NULL
      if(is.null(Report_data$d) | input$choice_visu == "dat"){
        File <- input$ydata_visu
        if(is.null(File))
          return(NULL)

        df <- rio::import(File$datapath)
      }
      else{
        if(input$bdata_visu == "Protein group"){
          df <- pg_ev$x
        }
        else if(input$bdata_visu == "Unique genes"){
          df <- gg_ev$x
        }
        else if(input$bdata_visu == "Peptides"){
          df <- pep_ev$x
        }
        else if(input$bdata_visu == "Peptides.MaxLFQ"){
          df <- peplfq_ev$x
        }
        else if(input$bdata_visu == "Precursors"){
          df <- precu_ev$x
        }
      }
      if(!is.null(df)){
        df <- as.data.frame(df)
        names(df)[1] <- "id"
      }
      df
    })
    output$visudata_up <- reactive({
      return(!is.null(visu_data()))
    })
    outputOptions(output, "visudata_up", suspendWhenHidden = FALSE)

    observe({
      if(!is.null(Report_data$d) & input$choice_visu == "base"){
        updateTextInput(session, "nmid_visu", value = "")
      }
    })
    output$top3_or_ibaq <- reactive({
      l <- FALSE
      dt <- "intensity"
      if(!is.null(visu_data())){
        l <- str_extract(colnames(visu_data()), "^Top3|^iBAQ")
        l <- l[!is.na(l)]
        dt <- c(dt, unique(l))
        l <- length(l) > 0
      }
      if(l){
        updateSelectInput(session, "dtype_visu", choices = c(dt, "all"))
      }
      return(l)
    })
    outputOptions(output, "top3_or_ibaq", suspendWhenHidden = FALSE)

    observe({
      if(!is.null(visu_data())){
        df <- visu_data()
        to_rm <- str_which(colnames(df), "nbTrypticPeptides|peptides_counts_all|^pep_count_")
        df <- df[,-to_rm]
        if(input$dtype_visu == "Top3"){
          df <- df[,str_which(colnames(df), "^Top3_")]
        }
        else if(input$dtype_visu == "iBAQ"){
          df <- df[,str_which(colnames(df), "^iBAQ_")]
        }
        else if(input$dtype_visu == "intensity"){
          idx <- str_which(colnames(df), "^iBAQ_|^Top3_")
          if(length(idx) > 0){
            df <- df[,-idx]
          }
        }
        n <- lapply(df, class)
        n <- sum(n == "numeric")
        updateSliderInput(session, "maxna", max = n)
      }
    })

    heat_ev <- reactiveValues(
      h = NULL
    )
    heat <- reactive({
      nm <- input$nmid_visu
      if(str_length(nm) == 0){
        nm <- NULL
      }
      heatmapDIA(visu_data(), input$transfo_visu, input$maxna, input$prval_visu, nm, data_type = input$dtype_visu)
    })
    observeEvent(input$seeheat_visu, {
      if(!is.null(visu_data())){
        if(str_length(input$nmid_visu) != 0){
          idx <- str_which(names(visu_data()), paste0("^", input$nmid_visu, "$"))
          if(purrr::is_empty(idx)){
            showNotification(paste("Please provide a valid column name. You enter :", input$nmid_visu,
                                   "and the column names are :", paste(names(visu_data()), collapse = ", "), "."),
                             type = "error", duration = 6)
          }
          else{
            showNotification("Get interactive heatmap", type = "message", duration = 4)
            heat_ev$h <- heat()
          }
        }
        else if(is.null(Report_data$d) | input$choice_visu == "dat"){
          showNotification("Don't forget to type a column name !", type = "error", duration = 5)
        }
        else{
          showNotification("Get interactive heatmap", type = "message", duration = 4)
          heat_ev$h <- heat()
        }
      }
      else{
        showNotification("Your data are NULL ! Start the calculation for the data you selected
                         or import a file", type = "error")
      }
    })
    output$heat_visu <- renderPlotly({
      heat_ev$h
    })

    ## DENSITY
    dens_ev <- reactiveValues(
      d = NULL
    )
    dens <- reactive({
      densityDIA(visu_data(), input$transfoD_visu, input$area_visu, input$titD_visu, data_type = input$dtype_visu)
    })
    observeEvent(input$seedens_visu, {
      if(!is.null(visu_data())){
        showNotification("Get density plot", type = "message", duration = 4)
        dens_ev$d <- dens()
      }
      else{
        showNotification("Your data are NULL ! Start the calculation for the data you selected
                         or import a file", type = "error")
      }
    })
    output$dens_visu <- renderPlot({
      dens_ev$d
    })
    output$down_dens <- downloadHandler(
      filename = function() {
        paste0("DensityPlot_dia_", Sys.Date(), ".png")
      },
      content = function(file){
        ggsave(file, plot = dens_ev$d, device = "png",
               width = 10, height = 7)
      }
    )

    ## MDS
    mds_ev <- reactiveValues(
      m = NULL
    )
    mds <- reactive({
      MDS_DIA(visu_data(), input$transfoM_visu, input$titM_visu, data_type = input$dtype_visu)
    })
    observeEvent(input$seemds_visu, {
      if(!is.null(visu_data())){
        cl <- lapply(visu_data(), class)
        cl <- cl == "numeric"
        if(sum(cl) >= 3){
          showNotification("Get MDS plot", type = "message", duration = 4)
          mds_ev$m <- mds()
        }
        else{
          showNotification("You need at list 3 columns in your data !", type = "error")
        }
      }
      else{
        showNotification("Your data are NULL ! Start the calculation for the data you selected
                         or import a file", type = "error")
      }
    })
    output$mds_visu <- renderPlot({
      mds_ev$m
    })
    output$down_mds <- downloadHandler(
      filename = function() {
        paste0("MDSPlot_dia_", Sys.Date(), ".png")
      },
      content = function(file){
        ggsave(file, plot = mds_ev$m, device = "png",
               width = 8, height = 8)
      }
    )


    ## PVV
    observe({
      if(stringr::str_length(input$designPVV_visu) != 0){
        updateSelectInput(session, "checkPVV_visu", choices = stringr::str_split(input$designPVV_visu, ",")[[1]])
      }
    })

    pvv_ev <- reactiveValues(
      m = NULL
    )
    pvv <- reactive({
      if(stringr::str_length(input$designPVV_visu) != 0){
        design <- stringr::str_split(input$designPVV_visu, ",")[[1]]
      }
      else{
        design <- NULL
      }
      validDIA(visu_data(), transformation = input$transfoPVV_visu, input$titPVV_visu, data_type = input$dtype_visu,
               design = design,
               to_check = input$checkPVV_visu, prop_cut = input$propcutPVV_visu)
    })
    observeEvent(input$seePVV_visu, {
      if(!is.null(visu_data())){
        showNotification("Get PVV plot", type = "message", duration = 4)
        pvv_ev$m <- pvv()
      }
      else{
        showNotification("Your data are NULL ! Start the calculation for the data you selected
                         or import a file", type = "error")
      }
    })
    output$PVV_visu <- renderPlot({
      pvv_ev$m
    })
    output$down_PVV <- downloadHandler(
      filename = function() {
        paste0("PVVPlot_dia_", Sys.Date(), ".png")
      },
      content = function(file){
        ggsave(file, plot = pvv_ev$m, device = "png",
               width = 8, height = 8)
      }
    )

    ## RT
    rt_ev <- reactiveValues(
      g = NULL,
      f = NULL
    )
    rtg <- reactive({
      g <- ggplot(Report_data$d, aes(iRT, RT, label1 = Precursor.Id, label2 = Protein.Ids, label3 = Genes, color = PG.Q.Value)) +
        geom_point() + facet_wrap(~File.Name) + labs(title = "Report data") + theme(plot.title = element_text(hjust = 0.5))
      ggplotly(g)
    })
    rtf <- reactive({
      d <- Report_data$d
      nm <- ""
      if(input$bdata_visu == "Protein group"){
        nm <- "Protein.Group"
      }
      else if(input$bdata_visu == "Unique genes"){
        nm <- "Genes"
      }
      else if(input$bdata_visu == "Peptides" | input$bdata_visu == "Peptides.MaxLFQ"){
        nm <- "Stripped.Sequence"
      }
      else if(input$bdata_visu == "Precursors"){
        nm <- "Precursor.Id"
      }
      d <- d[d[[nm]] %in% visu_data()$id,]
      g <- ggplot(d, aes(iRT, RT, label1 = Precursor.Id, label2 = Protein.Ids, label3 = Genes, color = PG.Q.Value)) +
        geom_point() + facet_wrap(~File.Name) + labs(title = "Report data filtered") + theme(plot.title = element_text(hjust = 0.5))
      ggplotly(g)
    })
    observeEvent(input$seert_visu, {
       showNotification("Get rentention time plot", type = "message", duration = 4)
       rt_ev$g <- rtg()
       rt_ev$f <- rtf()
    })
    output$rt1_visu <- renderPlotly({
      rt_ev$g
    })
    output$rt2_visu <- renderPlotly({
      rt_ev$f
    })

    ## PROTEOTYPIC
    ptyp_ev <- reactiveValues(
      g = NULL,
      f = NULL
    )
    ptypg <- reactive({
      ptyp <- Report_data$d[, c("File.Name", "Proteotypic")]
      ptyp$Proteotypic <- as.character(ptyp$Proteotypic)
      ggplot(ptyp, aes(Proteotypic, fill = Proteotypic)) +
        geom_bar() +
        facet_wrap(~File.Name) +
        labs(title = "Report data") +
        theme(plot.title = element_text(hjust = 0.5))
    })
    ptypf <- reactive({
      d <- Report_data$d
      nm <- ""
      if(input$bdata_visu == "Protein group"){
        nm <- "Protein.Group"
      }
      else if(input$bdata_visu == "Unique genes"){
        nm <- "Genes"
      }
      else if(input$bdata_visu == "Peptides" | input$bdata_visu == "Peptides.MaxLFQ"){
        nm <- "Stripped.Sequence"
      }
      else if(input$bdata_visu == "Precursors"){
        nm <- "Precursor.Id"
      }
      ptyp <- d[d[[nm]] %in% visu_data()$id,c("File.Name", "Proteotypic")]
      ptyp$Proteotypic <- as.character(ptyp$Proteotypic)

      ggplot(ptyp, aes(Proteotypic, fill = Proteotypic)) +
        geom_bar() +
        facet_wrap(~File.Name) +
        labs(title = "Report data filtered") +
        theme(plot.title = element_text(hjust = 0.5))
    })
    observeEvent(input$seeptyp_visu, {
      showNotification("Get proteotypic proportion", type = "message", duration = 4)
      ptyp_ev$g <- ptypg()
      ptyp_ev$f <- ptypf()
    })
    output$ptyp1_visu <- renderPlot({
      ptyp_ev$g
    })
    output$ptyp2_visu <- renderPlot({
      ptyp_ev$f
    })


    ### SELECT BEST WINDOWS
    shinyFileChoose(input, "replib_wsel", roots = volumes, session = session)
    output$filepaths2 <- renderPrint({
      if (is.integer(input$replib_wsel)) {
        cat("No files have been selected (shinyFileChoose)")
      } else {
        parseFilePaths(volumes, input$replib_wsel)
      }
    })

    Lib_data <- reactive({
      if (is.integer(input$replib_wsel)) {
        return(NULL)
      }
      else {
        File <- parseFilePaths(volumes, input$replib_wsel)
      }
      File <- parseFilePaths(volumes, input$replib_wsel)
      if(is.null(File))
        return(NULL)

      showNotification("Getting your data, this may take a while.", type = "message")
      diann_load(File$datapath)
    })
    lib_data <- reactiveValues(
      d = NULL
    )
    observe({
      lib_data$d <- Lib_data()
    })
    output$libdata_up <- reactive({
      return(!is.null(lib_data$d))
    })
    outputOptions(output, "libdata_up", suspendWhenHidden = FALSE)

    wsel_ev <- reactiveValues(
      o_h = NULL,
      o_hf = NULL,
      n_h = NULL,
      n_hf = NULL,
      both = NULL,
      both_f = NULL
    )
    wsel_result <- reactive({
      withCallingHandlers({
        shinyjs::html("bstw_wsel", "")
        get_bestwind(lib_data$d, bins = input$bins_wsel, per_frac = input$perfrac_wsel)
        },
      message = function(m) {
        shinyjs::html(id = "bstw_wsel", html = paste(m$message, "<br>", sep = ""), add = FALSE)
      }
      )
    })

    observeEvent(input$go_wsel, {
      showNotification("Get best windows for your data", type = "message", duration = 4)
      res <- wsel_result()
      wsel_ev$o_h <- res$orig_hist
      wsel_ev$o_hf <- res$orig_hist_perfrac
      wsel_ev$n_h <- res$new_hist
      wsel_ev$n_hf <- res$new_hist_perfrac
      wsel_ev$both <- ggpubr::ggarrange(res$orig_hist, res$new_hist, ncol = 2, nrow = 1)
      wsel_ev$both_f <- ggpubr::ggarrange(res$orig_hist_perfrac, res$new_hist_perfrac, ncol = 2, nrow = 1)
    })
    output$hist_wsel <- renderPlot({
      wsel_ev$both
    })
    output$hist_wsel_frac <- renderPlot({
      wsel_ev$both_f
    })
    output$downhist_wsel <- downloadHandler(
      filename = function() {
        paste0("DistributionWindow_dia_", Sys.Date(), ".png")
      },
      content = function(file){
        ggsave(file, plot = wsel_ev$both, device = "png",
               width = 12, height = 8)
      }
    )
    output$downhist_wsel_frac <- downloadHandler(
      filename = function() {
        paste0("DistributionWindowFraction_dia_", Sys.Date(), ".png")
      },
      content = function(file){
        ggsave(file, plot = wsel_ev$both_f, device = "png",
               width = 12, height = 8)
      }
    )
}


shinyApp(ui, server)







