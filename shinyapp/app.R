#############
#### App ####
#############

###############################################################################
#### Libraries and Sources ####################################################
###############################################################################
source("global.R")
# TODO: disable all warnings 
options(warn=-1)

#########################
#### Help Structures ####
#########################

# Table 
table_custom <- function(df, title = "", class = "", style = "", align = NULL, 
                         show_rownames_val = TRUE, show_colnames_val = TRUE) {
  
  if(show_colnames_val){
    header <- colnames(df)
  } else {
    header <- c()
  }
  
  if(show_rownames_val){
    show_rownames <- length(rownames(df)) > 0
  } else {
    show_rownames <- FALSE
  }
  
  
  if (is.null(align)) {
    align <- rep("l", ifelse(show_rownames, ncol(df) + 1, ncol(df)))
  }
  align[align == "l"] <- "text-start"
  align[align == "r"] <- "text-end"
  align[align == "c"] <- "text-center"
  
  if (show_rownames) {
    df$.rownames <- rownames(df)
  }
  div(
    class = class,
    style = style,
    HTML('<table class="table table-striped table-bordered"><thead><tr>'),
    HTML(ifelse(show_rownames, paste0('<th scope="col" class="text-primary ',align[1], '">', title,'</th>'), '')),
    HTML(paste0('<th scope="col" class="', ifelse(show_rownames, align[2:length(align)], align),'">',header , '</th>', collapse = '')),
    HTML('</tr></thead><tbody>'),
    HTML(paste0(
      '<tr>',
      apply(df, 1, function(row) {
        if (show_rownames) {
          paste0(
            '<th class="', align[1], '" scope="row">', row['.rownames'], '</th>',
            paste0('<td class="', align[2:length(align)], '">', row[which(names(row) != '.rownames')], '</td>', collapse = '')
          )
        } else {
          paste0(
            paste0('<td class="', align, '">', row, '</td>', collapse = '')
          )
        }
      }),
      '</tr>',
      collapse = ''
    )),
    HTML('</tbody></table>')
  )
}

# Fancy table 
fancyDT <- function(dt,filename,...){
  tmp_tbl <- datatable(dt, 
            # style = 'bootstrap',
            # class = 'cell-border',
            escape = FALSE, 
            filter = 'top', 
            extensions = c('KeyTable'),
            selection = 'none',
            
            options = list(
              dom = '<"top"B>frit<"bottom"lp>',
              paging = TRUE,
              # pageLength = 5, 
              lengthMenu = list(c(5,10,25,50, -1), c("5","10","25", "50","All")), 
              scroller = TRUE,
              scrollX = TRUE,
              autoWidth = FALSE,
              keys = TRUE,
              select = TRUE,
              
              
              rownames= FALSE,
              # server = FALSE,
              ... # additional options passed with function
              
              
            ))
  
}




#####################
#### UI elements ####
#####################
ui <- dashboardPage(
  
  
  # Dashboard configuration
  skin =  "red", 
  title = "NeEDL",
  freshTheme = NULL,
  preloader = NULL,
  md = FALSE,
  scrollToTop = FALSE,
  options = list(sidebarExpandOnHover = TRUE),
  
  
  
  # Header
  header = dashboardHeader(
    
    title = tagList(
      span(class = "logo-lg", div(img(src = needl_img_path, style= "width:20%"), "NeEDL")), 
      img(src = needl_img_path, style= "width:150%")
      ),
    
    
    titleWidth = NULL,
    disable = FALSE,
    .list = NULL,
    leftUi = NULL,
    controlbarIcon = shiny::icon("question"),
    fixed = FALSE
    # dropdownMenu(type = c("messages", "notifications", "tasks"),
    #              badgeStatus = "primary", icon = NULL, headerText = NULL, .list = NULL)
  ),
  
  
  
  # Sidebar left
  sidebar = dashboardSidebar(
    id = "sidebar",
    disable = FALSE,
    width = NULL,
    collapsed = FALSE,
    minified = TRUE,
    
    
    sidebarMenu(
      style = "position:fixed;width: inherit;",
      id = NULL, 
      .list = NULL, 
      
      menuItem(text = "General", 
               icon = icon("home"), 
               badgeLabel = NULL, 
               badgeColor = "green",
               tabName = "general", 
               href = NULL, 
               newtab = TRUE, 
               selected = NULL,
               startExpanded = FALSE
               ), 
      
      menuItem(text = "Technical details", 
               icon = icon("link"), 
               badgeLabel = NULL, 
               badgeColor = "green",
               tabName = "technical", 
               href = NULL, 
               newtab = TRUE, 
               selected = NULL,
               startExpanded = FALSE
               ), 
      
      menuItem(text = "Epistasis disease atlas", 
               icon = icon("book-atlas"), 
               badgeLabel = NULL, 
               badgeColor = "green",
               tabName = "epiatlas", 
               href = NULL, 
               newtab = TRUE, 
               selected = NULL,
               startExpanded = FALSE
               )
      )
    ),
  
  
  # Body
  body = dashboardBody(
    # For nested box content 
    tags$head(tags$script(HTML("$(() => $('body').on('shown.bs.collapse', '.box', function(evt) { 
      setTimeout(function(){
         $(evt.target).find('.shiny-bound-output').trigger('shown.bs.collapse');
      }, 800);
   }))"))),
   
   tags$script("$(document).on('shiny:connected', function(event) {
                  var myWidth = $(window).width();
                  Shiny.onInputChange('shiny_width',myWidth)
                  
                  });"),
          
    tags$script("$(document).on('shiny:connected', function(event) {
                  var myHeight = $(window).height();
                  Shiny.onInputChange('shiny_height',myHeight)
                  
                  });"),

   shinyjs::useShinyjs(),
   shinyFeedback::useShinyFeedback(),
   
    tabItems(
      ##################
      ##### General ####
      ##################
      tabItem(tabName = "general", 
              
              # Heading
              h2(
                HTML('<i class="fa-solid fa-home"></i>'),
                "General"
              ),
              
              
              box(
                width = 12, 
                title = "About NeEDL", 
                solidHeader = TRUE, 
                # background = "red", 
                status = "danger", 
                collapsible = TRUE,
                collapsed = TRUE,
              
    
                h4(align = "justify",
                  style = "font-size: 500; line-height: 1.5;",
                  HTML(
                    "Genome-wide association studies (GWAS) aim to link genetic variants to phenotypic traits of interest, most commonly a disease. More specifically, GWAS usually looks for biallelic single nucleotide polymorphisms (SNPs) that are individually predictive of the phenotype. While thousands of individual SNPs have been associated with diseases since the early 2000s, they account only for a fraction of the investigated traits’ heritability. The most common hypothesis is that the missing heritability can be explained by epistasis, i.e., by interactions between SNPs that are jointly predictive of the phenotype but individually have little or no effect. However, although epistasis is assumed to play an important role in the genomics of complex phenotypic traits, no undisputed cases of epistasis in humans are known. Therefore, scalable epistasis detection tools that yield interpretable, high-quality results are urgently needed to gain further knowledge of possible genetic causes of diseases.", 
                    "<br> Developing epistasis detection tools that can capture biologically interpretable candidate SNP sets is difficult for at least four reasons: firstly, the problem of epistasis detection needs to be modeled formally in order to render it algorithmically accessible, and there seems to be no consensus in the biomedical community as to the choice of the formal model. Secondly, it is often unclear whether predicted cases of epistasis are biologically meaningful or mere statistical artifacts. Thirdly, especially if interactions of order higher than two are to be considered, the increasing availability of large-scale GWAS data makes epistasis detection computationally expensive. Existing frequently used epistasis detection tools such as Potpourri, LinDen, PoCos, MACOED, and BiologicalEpistasis therefore do not scale to large datasets and are mostly restricted to pairwise interactions. Fourthly, not everyone has the computational power to carry out such expensive computations, and there is a lack of easily accessible resources that allow medical researchers to browse and explore pre-computed epistasis candidates in an interactive way.", 
                    "<br> To address these challenges, we here present NeEDL, a quantum-computing-ready network-based epistasis detection tool based on local search. Thanks to its modular and easy-to-use design, NeEDL can be configured to use various statistical epistasis models10 that associate higher-order genotypes with phenotype data, which addresses the first of the four above-mentioned challenges. To address the second challenge and prioritize biologically plausible SNP sets, NeEDL makes use of an SNP-SNP interaction (SSI) network, which we construct by mapping SNPs to proteins (using dbSNP) and then connecting two SNPs if they affect the same or neighboring proteins in a protein-protein interaction (PPI) network. Here, our hypothesis is that when multiple SNPs all affect the same protein or multiple interacting proteins responsible for a specific function, they are more likely to be involved in epistatic interaction, as it might happen that multiple, individually innocuous changes to a protein (complex) lead to a loss-of-function scenario. Within the SSI network, NeEDL uses local search with multi-start and simulated annealing to find connected subgraphs of a user-specified size that are locally optimal w. r. t. the selected statistical epistasis model. Focusing on SNP sets inducing connected subgraphs in the SSI network not only increases the likelihood of uncovering biologically meaningful cases of epistasis but also dramatically reduces the size of the search space. To further decrease the runtime and ensure scalability to large-scale GWAS data, we provide NeEDL in a highly parallelized C++ implementation. Finally, NeEDL incorporates quantum computing (QC) algorithms that allow computing high-quality initial solutions for the local search and will hence further reduce NeEDL’s runtime requirements should production-ready QC hardware become available in the future."
                    )
                  )
                ),
              
              # Workflow of NeEDL
              box(
                width = 12, 
                title = "Workflow of NeEDL", 
                solidHeader = TRUE, 
                # background = "red", 
                status = "danger",
                collapsible = TRUE,
                collapsed = TRUE,
                
                # Image of NeEDL workflow
                div(
                  img(src = needl_workflow_img_path, style = 'width:100%'),
                  style= 'max-width:800px;margin:0 auto'
                )
              )
              ), 
      
      ############################
      ##### Technical Details ####
      ############################
      tabItem(tabName = "technical", 
              
              # Heading
              h2(
                HTML('<i class="fa-solid fa-link"></i>'),
                "Technical details"
              ),
              
              # Parameters of the run
              box(
                width = 12, 
                title = "Pipeline configuration",  
                solidHeader = TRUE, 
                # background = "red", 
                status = "danger",
                collapsible = TRUE,
                collapsed = TRUE,
                
                h4(align = "justify",
                  style = "font-size: 500; line-height: 1.5;",
                  HTML(
                   "The following table offers a comprehensive view of each stage
                   of the NeEDL pipeline, including the corresponding variables 
                   and their respective values for each step."
                  )
                ),
                
                # Download button for pipeline configuration
                downloadBttn("pipeline_download", 
                             "Download", 
                             style = "unite",
                             color = "primary",
                             size = "sm",
                             block = FALSE,
                             no_outline = TRUE),
              
                h5(' '), 
                
                
                # Table with pipeline configuration
                div(htmlOutput("task_stats_tbl"))
               
                
              ),
              
              # Runtime information
              fluidRow(
                column(5,
                  box(
                    width = 12, 
                    title = "Runtime",  
                    solidHeader = TRUE, 
                    # background = "red", 
                    status = "danger",
                    collapsible = TRUE,
                    collapsed = TRUE,
                
                    h4(align = "justify",
                       style = "font-size: 500; line-height: 1.5;",
                       HTML(
                         "The following table provides runtime information 
                         of the NeEDL run. "
                       )
                    ),
                   
                    # Download button for runtime information
                    downloadBttn("runtime_download", 
                                 "Download", 
                                 style = "unite",
                                 color = "primary",
                                 size = "sm",
                                 block = FALSE,
                                 no_outline = TRUE),
                    
                    h5(' '), 
                    
                    # time parameters of the run
                    div(htmlOutput("task_time_tbl"))
                  )
                ),
                
                column(7, 
                  # Memory information
                  box(
                    width = 12, 
                    title = "Memory",  
                    solidHeader = TRUE, 
                    # background = "red"
                    status = "danger",
                    collapsible = TRUE,
                    collapsed = TRUE,
                    
                    h4(align = "justify",
                       style = "font-size: 500; line-height: 1.5;",
                       HTML(
                         "In the following table the memory requirements of 
                         the NeEDL run are presented. "
                       )
                    ),
                    
                    # Download button for memory information
                    downloadBttn("memory_download", 
                                 "Download", 
                                 style = "unite",
                                 color = "primary",
                                 size = "sm",
                                 block = FALSE,
                                 no_outline = TRUE),
                    
                    h5(' '), 
                    
                    # memory parameters of the run
                    div(htmlOutput("task_memory_tbl"))
                  )
                )
              )
      
      ), 
      
      
      ###################################################################
      ##### Epistasis disease atlas with network centrality analysis #### 
      ###################################################################
      tabItem(tabName = "epiatlas", 
              
              
              # Heading 
              h2(
                img(src = epiatlas_img_path, style= "width:6%"), 
                "Epistasis disease atlas"
              ),
              
              # InfoBox for information on whole network and seeding
              infoBox(width = 12,
                      title = h4("Information on whole SSI network and seeding"), 
                      subtitle = h5("The following section offers general information
                                    about the entire SNP-SNP interaction network."), 
                      icon = icon("1"),
                      color = "navy", 
                      fill = TRUE
              ),
              
              
              # Network statistics
              box(
                width = 12, 
                title = "Network statistics",
                solidHeader = TRUE, 
                # background = "red"
                status = "primary",
                collapsible = TRUE,
                collapsed = TRUE,
                
                h4(align = "justify",
                   style = "font-size: 500; line-height: 1.5;",
                   HTML(
                     'The following table displays information regarding the 
                     number of nodes and edges present in the network. 
                     The edges are classified based on their origin and 
                     corresponding numbers are provided. 
                     The origin pertains to the underlying foundation for 
                     linking two SNPs during network construction.
                     If these SNPs affect the same proteins, the edge connecting
                     them is labeled as "SAME TAG". Conversely, if they affect 
                     neighboring proteins in the PPI network, the edge is 
                     labeled as "BIOGRID".'
                   )
                ),
                
                # Download button for network statistics
                downloadBttn("net_stats_download", 
                             "Download", 
                             style = "unite",
                             color = "primary",
                             size = "sm",
                             block = FALSE,
                             no_outline = TRUE),
                h5(' '), 
                
                
                # Network statistics 
                div(shinycssloaders::withSpinner(
                  htmlOutput("net_stats_tbl_epi")
                  )
                )
              ), 
              
              box(
                width = 12, 
                title = "Distributions",
                solidHeader = TRUE, 
                # background = "red"
                status = "primary",
                collapsible = TRUE,
                collapsed = TRUE,
                
                # Score distributions 
                tags$style('.nav-tabs-custom .nav-tabs li.active {
                          border-top-color: #d73925;
                      }"',
                ),
                
                h4(align = "justify",
                   style = "font-size: 500; line-height: 1.5;",
                   HTML(
                     "The distributions of four network centralities across the
                     whole SSI network are provided, including degree, eigenvector, 
                     page rank and closeness. However, due to the increased 
                     computational complexity when calculating all distances 
                     to other nodes for each node, the closeness distribution 
                     is approximated based on one thousand randomly drawn data 
                     points from both candidate and non-candidate SNPs. 
                     The distributions are displayed as histograms with the x-axis 
                     representing the corresponding centrality, and the y-axis 
                     showing the logarithmically transformed counts to base ten."
                   )
                ),
                
                fluidRow(
                  column(3, 
                         # Distribution of selected centrality in the SSI network
                         selectInput("dist_select", 
                                     label = h5("Select the centrality to display"), 
                                     choices = list("Degree" = "degree",
                                                    "Eigenvector" = "eigenvector",
                                                    "Pagerank" = "pagerank", 
                                                    "Closeness" = "closeness"
                                     ), 
                                     selected = "degree"
                         )
                  )
                ),
                
                
                h5(' '), 
                
                # Download button for centrality distribution data
                downloadBttn("centr_dist_data_download", 
                             "Download data", 
                             style = "unite",
                             color = "primary",
                             size = "sm",
                             block = FALSE,
                             no_outline = TRUE),
                
                h5(' '), 
                
                # Download button for centrality distribution plots
                downloadBttn("centr_dist_plt_download", 
                             "Download plot", 
                             style = "unite",
                             color = "primary",
                             size = "sm",
                             block = FALSE,
                             no_outline = TRUE),
                
                h5(' '), 
             
                shinycssloaders::withSpinner(
                  plotOutput(outputId = "dist_centr_plt", height = "500px")
                )
                
              ), 
             
              # Start seeds 
              box(
                width = 12, 
                title = "Start seeds",
                solidHeader = TRUE, 
                # background = "red"
                status = "primary",
                collapsible = TRUE,
                collapsed = TRUE,
                
                h4(align = "justify",
                   style = "font-size: 500; line-height: 1.5;",
                   HTML(
                     "The following searchable and sortable table provides 
                     details on the initial seed pairs (RS_IDS). 
                     For each seed pair the table includes 14 different scores
                     based on different statistical models. The table is ranked 
                     based on one of these scores. Additionally, the origin of 
                     the seeds is displayed, indicating which seeding routine 
                     was used to select the corresponding start seeds."
                   )
                ),
                
                # Download button for network statistics
                downloadBttn("seeds_tbl_download", 
                             "Download", 
                             style = "unite",
                             color = "primary",
                             size = "sm",
                             block = FALSE,
                             no_outline = TRUE),
                
                
                rclipboardSetup(),
                # Table with information on start seeds
                shinycssloaders::withSpinner(
                  DTOutput("seeds_tbl_epi")
                )
              ), 
              
              
              # Info box for analysis of candidate SNP sets
              infoBox(width = 12,
                      title = h4("Analysis of candidate SNP sets"), 
                      subtitle = h5("In this section the candidate SNP sets can 
                                    be analyzed using the set scores, the SNP-SNP 
                                    interaction network and several network centralities."), 
                      icon = icon("2"),
                      color = "navy", 
                      fill = TRUE
              ),
              
              
              # Search score over time
              box(
                width = 12, 
                title = "Search score over time",
                solidHeader = TRUE, 
                # background = "red"
                status = "primary",
                collapsible = TRUE,
                collapsed = TRUE,
                
                
                h4(align = "justify",
                   style = "font-size: 500; line-height: 1.5;",
                   HTML(
                     "To examine the score development of the candidate SNP sets, 
                     we provide a step plot illustrating the trajectory of the 
                     Maximum Likelihood Model (MLM) score over time in milliseconds."
                   )
                ),
                
                
                # Download button for search score over time data
                downloadBttn("search_score_data_download", 
                             "Download data", 
                             style = "unite",
                             color = "primary",
                             size = "sm",
                             block = FALSE,
                             no_outline = TRUE),
                
                h5(' '), 
                
                # Download button for search score over time plot
                downloadBttn("search_score_plt_download", 
                             "Download plot", 
                             style = "unite",
                             color = "primary",
                             size = "sm",
                             block = FALSE,
                             no_outline = TRUE),
                
                h5(' '), 
                
                # Search score over time plot
                shinycssloaders::withSpinner(
                  plotOutput(outputId = "search_score_plt", height = "500px")
                ) 
                
              ), 
              
              # Set scores 
              box(
                width = 12,  
                title = "Set scores",
                solidHeader = TRUE, 
                # background = "red"
                status = "primary",
                collapsible = TRUE,
                collapsed = TRUE,
                
                
                h4(align = "justify",
                   style = "font-size: 500; line-height: 1.5;",
                   HTML(
                     "To further investigate the distribution of the 14 previous 
                     presented scores, we offer a boxplot displaying the 
                     distribution of each score across all candidate SNP sets."
                   )
                ),
                
                
                # Download button for search score over time plot
                downloadBttn("set_scores_data_download", 
                             "Download data", 
                             style = "unite",
                             color = "primary",
                             size = "sm",
                             block = FALSE,
                             no_outline = TRUE),
                
                h5(' '), 
                
                # Download button for set scores
                downloadBttn("set_scores_plt_download", 
                             "Download plot", 
                             style = "unite",
                             color = "primary",
                             size = "sm",
                             block = FALSE,
                             no_outline = TRUE),
                
                # Set scores 
                shinycssloaders::withSpinner(
                  plotOutput(outputId = "set_scores_plt", height = "500px")
                )
              ),
              
              # Candidate Sets
              box(
                width = 12, 
                title = "Candidate sets",
                solidHeader = TRUE, 
                # background = "red"
                status = "primary",
                collapsible = TRUE,
                collapsed = TRUE,
                
                h4(align = "justify",
                   style = "font-size: 500; line-height: 1.5;",
                   HTML(
                     'Analogously to the searchable and sortable start seed table, 
                     a table for the candidate SNP sets (RS_IDS) is available, 
                     except for the origin column. 
                     Additionally, the table includes a "copy" button located in
                     the right column of the rsIDs of the candidate SNP sets, 
                     which allows for copying the rsIDs in the correct format 
                     for subsequent analysis.'
                   )
                ),
                
                # Download button for candidate table
                downloadBttn("cand_tbl_download", 
                             "Download", 
                             style = "unite",
                             color = "primary",
                             size = "sm",
                             block = FALSE,
                             no_outline = TRUE),
                
                rclipboardSetup(),
                shinycssloaders::withSpinner(
                  DTOutput("cand_tbl_epi")
                )
              ), 
              
              # SNP-SNP interaction network
              box(
                width = 12, 
                title = "SNP-SNP interaction network",
                solidHeader = TRUE, 
                # background = "red"
                status = "primary",
                collapsible = TRUE,
                collapsed = TRUE,
                
                h4(align = "justify",
                   style = "font-size: 500; line-height: 1.5;",
                   HTML(
                     'To begin exploring the disease data, a set of SNPs can be 
                     entered into the query field (separated by semicolons, commas, or whitespaces).<br>
                     For the construction of the corresponding SSI network, 
                     the immediate neighborhood of the selected SNPs can be 
                     filtered using network centralities, including degree, 
                     eigenvector, and page rank. For each of them a range can 
                     be selected as a prerequisite for the centralities of the
                     neighboring SNPs. The selected centralities can be prioritized 
                     by drag and drop, with the most important centrality at the top. 
                     However, due to the computational burden mentioned before, 
                     the closeness centrality is not provided.<br> 
                     After selecting at least one centrality, it is possible 
                     to choose the number of neighboring SNPs to display that
                     meet the parameter combinations. After requesting the 
                     network construction via the “Create Network” button, 
                     the network is built, and results can be further explored.'
                   )
                ),
                
                # Box for input
                box(
                  width = 12, 
                  title = "Input",
                  solidHeader = FALSE,
                  status = "primary",
                  collapsible = FALSE,
                  
                  
                # rs-ids to visualize network
                textInput(inputId = "selected_SNPs_net_epi", 
                          # label = "Enter SNPs:", 
                          h5("Enter rsIDs for SNPs:"),
                          width = "100%",
                          placeholder = "rs7596121;rs41440544;rs41501252;rs16976638;rs16976644;rs10518828;rs16976648"
                          # "rs7596121;rs41440544;rs41501252;rs16976638;rs16976644;rs10518828;rs16976648"
                ),
                
                
                fluidRow(
                  column(8,
                    fluidRow(
                      # Centralities to filter for
                      column(4, 
                             # Degree
                             selectInput("centr_filter_degree", 
                                         label = h5("Degree"), 
                                         choices = list("ascending" = 1, 
                                                        "descending" = -1, 
                                                        "none" = 0), 
                                         selected = 0)
                             ),
                      column(8,                              
                             # Range for degree
                             sliderInput(inputId = "degree_slider", 
                                         label = h5("Select min and max degree:"),
                                         min = min(nodes$degree), 
                                         max = max(nodes$degree), 
                                         value = c(min(nodes$degree), max(nodes$degree))
                             )
                      )
                    ),
                    
                    
                    
                    fluidRow(
                      column(4, 
                             # Eigenvector
                             selectInput("centr_filter_eigenvector", 
                                         label = h5("Eigenvector"), 
                                         choices = list("ascending" = 1, 
                                                        "descending" = -1, 
                                                        "none" = 0), 
                                         selected = 0)
                      ), 
                      column(8,                              
                             # Range for eigenvector
                             sliderInput(inputId = "eigenvector_slider", 
                                         label = h5("Select min and max eigenvector:"),
                                         min = 0, 
                                         max = 1, 
                                         value = c(0, 1)
                                         # min = min(nodes$eigenvector), 
                                         # max = max(nodes$eigenvector), 
                                         # value = c(min(nodes$eigenvector), max(nodes$eigenvector))
                             )
                      )
                    ),
                    
                    fluidRow(
                      column(4, 
                             # Page rank
                             selectInput("centr_filter_pagerank", 
                                         label = h5("Page rank"), 
                                         choices = list("ascending" = 1, 
                                                        "descending" = -1, 
                                                        "none" = 0), 
                                         selected = 0)
                      ), 
                      column(8,                              
                             # Range for pagerank
                             sliderInput(inputId = "pagerank_slider", 
                                         label = h5("Select min and max page rank:"),
                                         min = 0, 
                                         max = 1, 
                                         value = c(0, 1)
                                         # min = min(nodes$pagerank), 
                                         # max = max(nodes$pagerank), 
                                         # value = c(min(nodes$pagerank), max(nodes$pagerank))
                             )
                      )
                    )
                  ),
                  
                  column(4, 
                         # Prioritization for network centralities 
                         uiOutput("centr_priority")
                  )
                ),
                
                
                
                fluidRow(
                  column(6, 
                         # number of candidates to show
                         sliderInput(inputId = "num_res", 
                                     label = h5("Number of non-candidate neighbours:"),
                                     min = 0, 
                                     max = 30, 
                                     value = 0
                         )
                  )
                ),
                
                fluidRow(
                  column(3, 
                         h1(" "),
                         # action button to start generating the network
                         shinyWidgets::actionBttn(inputId = "snp_net_action",
                                                    label = "Create Network", 
                                                    icon = NULL,
                                                    style = "unite",
                                                    color = "primary",
                                                    size = "sm",
                                                    block = FALSE,
                                                    no_outline = TRUE),
                         h1(" ")
                  ) 
                )
                ),
                
                h4(align = "justify",
                   style = "font-size: 500; line-height: 1.5;",
                   HTML(
                     'In the SSI network, the yellow stars represent the 
                     candidate SNPs, while the blue circles represent the 
                     neighbors that were chosen based on the selected centralities. 
                     If no neighbors are displayed or less than the specified number,
                     it means that either no neighbors were found or not enough 
                     neighbors met the specified conditions for the given parameter 
                     combination. In the network, the size of a node is 
                     proportional to its degree, i.e., a larger node indicates 
                     a higher degree. <br>
                     To choose a specific rsID, users can make use of the 
                     "Select by id" field located in the top left corner of 
                     the network. This will highlight the chosen rsID in the network.'
                   )
                ),
                
                fluidRow(
                  column(12, 
                         # network visualization
                         shinycssloaders::withSpinner(
                           visNetworkOutput(outputId = "network_epi",
                                            width = "100%", 
                                            height = "650px")
                           )
                         )
                ),
                
                
                h5(' '), 
                
                # Download button for current SSI network
                downloadBttn("snp_net_download", 
                             "Download network", 
                             style = "unite",
                             color = "primary",
                             size = "sm",
                             block = FALSE,
                             no_outline = TRUE),
                
                h5(' '),
                
                
                fluidRow(
                  column(8, 
                         # Table with information about the current SNP
                         box(
                           width = 12, 
                           title = "SNP information", 
                           solidHeader = TRUE, 
                           # background = "red"
                           status = "primary",
                           collapsible = TRUE,
                           collapsed = FALSE,
                           
                           h4(align = "justify",
                              style = "font-size: 500; line-height: 1.5;",
                              HTML(
                                "By clicking on an SNP in the network, corresponding
                                information about it can be obtained, such as its
                                rsID, its interacting gene or genes, and its 
                                network centralities. For SNPs not contained in 
                                any candidate SNP set the degree, eigenvector, 
                                and page rank are displayed. For candidate SNPs, 
                                the closeness centrality, the 14 statistical scores,
                                and the rank according to one of the scores are 
                                available as well.<br>
                                However, the time required for this process may
                                vary depending on the number of displayed SNPs, 
                                as it entails searching for the interacting genes 
                                in BIOGRID for each gene."
                              )
                           ),
                     
                           # Table with information on given SNP
                           shinycssloaders::withSpinner(
                             htmlOutput("snp_tbl_epi")
                           ) 
                           
                         ),
                  ), 
                  column(4, 
                         # Network statistics for current SNP-SNP interaction network
                         box(
                           width = 12, 
                           title = "Network statistics",
                           solidHeader = TRUE, 
                           # background = "red"
                           status = "primary",
                           collapsible = TRUE,
                           collapsed = FALSE,
                           
                           h4(align = "justify",
                              style = "font-size: 500; line-height: 1.5;",
                              HTML(
                                "Similar to the network statistics of the entire 
                                network, the following table provides the number 
                                of nodes and edges for the SSI network based on
                                the selected SNPs, together with the numbers for
                                the different edge types."
                              )
                           ),
                           
                           # Download button for network statistics
                           downloadBttn("snp_net_stats_download", 
                                        "Download", 
                                        style = "unite",
                                        color = "primary",
                                        size = "sm",
                                        block = FALSE,
                                        no_outline = TRUE),
                           h5(' '), 
                           
                           # Network statistics 
                             div(htmlOutput("snp_net_stats_tbl_epi"))
                             
                           ), 
                         # SNP databases
                         box(
                           width = 12, 
                           title = "SNP databases",
                           solidHeader = TRUE, 
                           # background = "red"
                           status = "primary",
                           collapsible = TRUE,
                           collapsed = FALSE,
                           
                           h4(align = "justify",
                              style = "font-size: 500; line-height: 1.5;",
                              HTML(
                                "The following three databases allow for the 
                                analysis of the currently viewed SNP. 
                                If no SNP is selected, the home page 
                                will be displayed."
                              )
                           ),
                           
                           fluidRow(
                             column(12, 
                                    # dbSNP
                                    uiOutput(outputId = "dbSNP_link_epi"),
                                    h4(' '),
                                    # ClinVar
                                    uiOutput(outputId = "clinvar_link_epi"), 
                                    h4(' '),
                                    # SNPmap
                                    uiOutput(outputId = "snpmap_link_epi")
                                    )
                             )
                           )
                         )
                )
              ),
              
              
              # Information on selected SNPs
              box(
                width = 12, 
                title = "Information on selected SNPs",
                solidHeader = TRUE, 
                # background = "red"
                status = "primary",
                collapsible = TRUE,
                collapsed = TRUE,
                
                uiOutput("info_on_selected_SNPs_box")
              ),
              
              
              # InfoBox for analysis of protein-protein interaction network
              infoBox(width = 12,
                      title = h4("Analysis of protein-protein interaction network"), 
                      subtitle = h5("In this section the genes/proteins corresponding to the 
                      SNPs in the SNP-SNP interaction network can be analyzed 
                      using different databases."), 
                      icon = icon("3"),
                      color = "navy", 
                      fill = TRUE
              ),
              # Protein-protein interaction network
              box(
                width = 12, 
                title = "Protein-protein interaction network",
                solidHeader = TRUE, 
                # background = "red"
                status = "primary",
                collapsible = TRUE,
                collapsed = TRUE,
                
                uiOutput("gene_gene_box")
              ), 
              
              # Genome browser box
              box(
                width = 12, 
                title = "Genome browser",
                solidHeader = TRUE, 
                # background = "red"
                status = "primary",
                collapsible = TRUE,
                collapsed = TRUE,
                
                # PlugIn for DrugStone using Genome Browser
                htmlOutput("genome_browser_frame")
              )
      )
    )
  ),
  
  # Controlbar 
  controlbar = dashboardControlbar(
    controlbarMenu(
      # h5('See Github for more information.')
      controlbarItem(
        title = "See Github for more information."
      )
      )
    ),
  
  
  

)



######################
#### Server logic ####
######################
server <- function(input, output, session) {
  
  
  ##################
  ##### General ####
  ##################
  ############################
  ##### Technical Details ####
  ############################
  
  # Parameters of the run from pipeline configuration
  output$task_stats_tbl <- renderText({
    as.character(
      div(
        table_custom(pipeline_config_dt,
                     title = "",
                     class="col-md-6",
                     show_rownames_val = FALSE,
                     show_colnames_val = FALSE),
      )
    )
  })

  # Pipeline download
  output$pipeline_download <- downloadHandler(
    filename <- "pipeline_config.csv",
    
    content = function(file) {
      write.csv(pipeline_config_dt, file, row.names = FALSE)
    }
  )
  
  
  # Time parameters of the run
  output$task_time_tbl <- renderText({
    as.character(
      div(
        table_custom(run_log_dt[1:3], 
                     title = "", 
                     class="col-md-12", 
                     show_rownames_val = FALSE),
      )
    )
  })
  
  # Runtime download
  output$runtime_download <- downloadHandler(
    filename <- "runtime.csv",
    
    content = function(file) {
      write.csv(run_log_dt[1:3], file, row.names = FALSE)
    }
  )
  
  # Time parameters of the run
  output$task_memory_tbl <- renderText({
    as.character(
      div(
        table_custom(run_log_dt[4:11], title = "", 
                     class="col-md-12", show_rownames_val = FALSE),
      )
    )
  })
  
  # Memory download
  output$memory_download <- downloadHandler(
    filename <- "memory.csv",
    
    content = function(file) {
      write.csv(run_log_dt[4:11], file, row.names = FALSE)
    }
  )
  
  
  ###################################################################
  ##### Epistasis disease atlas with network centrality analysis ####
  ###################################################################

  
  #############################################################
  # EDA: PART 1: INFORMATION ON WHOLE SSI NETWORK AND SEEDING #
  #############################################################
  
  ### Network stats ###
  output$net_stats_tbl_epi <- renderText({
    as.character(
      div(
        table_custom(get_net_stats(nodes, edges), title = "", class="col-md-5", 
                     show_rownames_val = FALSE)
      )
    )
  })

  # Network stats download
  output$net_stats_download <- downloadHandler(
    filename <- "network_stats.csv",
    
    content = function(file) {
      write.csv(get_net_stats(nodes, edges), file, row.names = FALSE)
    }
  )
  
  
  
  
  ### Distribution of centralities ###
 
  # Centrality distribution plot
  dist_centr_plt_re <- reactive({
    plot_centrality(nodes, centrality = input$dist_select, num_bins = 50)
  })
  
  
  # Distribution of selected centrality in the SSI network
  output$dist_centr_plt <- renderPlot({
    dist_centr_plt_re()
  })
  
  # Download for distribution of selected centrality data
  output$centr_dist_data_download <- downloadHandler(
    filename = function() { paste(input$dist_select, '_centr_dist.csv', sep='') },
    
    content = function(file) {
      tmp_nodes <- nodes[, .(rs_id, get(input$dist_select))]
      colnames(tmp_nodes) <- c('rsID', input$dist_select)
      write.csv(tmp_nodes, file, row.names = FALSE)
    }
  )
  

  # Download for distribution of selected centrality plot
  output$centr_dist_plt_download <- downloadHandler(
    filename = function() { paste(input$dist_select, '_centr_dist.png', sep='') },

    content = function(file) {
      png(file,
          width = input$shiny_width,
          height = input$shiny_height)

      print(dist_centr_plt_re())
      dev.off()
    }
  )
  
  
  
  ## Start seeds ##
  
  # Output for start seeds table 
  output$seeds_tbl_epi <- DT::renderDT({
    fancyDT(seeds_dt, "start_seeds")
  })
  
  # Network stats download
  output$seeds_tbl_download <- downloadHandler(
    filename <- "seeds.csv",
    
    content = function(file) {
      write.csv(seeds_dt, file, row.names = FALSE)
    }
  )
  
  
  ##############################################
  # EDA: PART 2: ANALYSIS OF CANDIDATE SNP SET #
  ##############################################
  
  ### Candidates ###
  
  ## Search score over time
  
  # Search score over time plot
  search_score_plt_re <- reactive({
    plot_search_score(search_score_dt)
  })
  
  # Search score over time 
  output$search_score_plt <- renderPlot({
    search_score_plt_re()
  })
  
  # Search score over time data download
  output$search_score_data_download <- downloadHandler(
    filename <- "search_score_over_time.csv",
    
    content = function(file) {
      write.csv(search_score_dt, file, row.names = FALSE)
    }
  )
  
  # Search score over time plot download
  output$search_score_plt_download <- downloadHandler(
    filename <- "search_score_over_time.png",
    
    content = function(file) {
      png(file,
          width = input$shiny_width,
          height = input$shiny_height)
      
      print(search_score_plt_re())
      dev.off()
    }
  )
  
  # Search score over time plot
  set_scores_plt_re <- reactive({
    plot_set_scores(res_dt = res_ext_dt, score_names = score_names)
  })
  
  ## Set scores 
  
  # Set scores plot
  output$set_scores_plt <- renderPlot({
    set_scores_plt_re()
  })
  
  
  # Set score data download
  output$set_scores_data_download <- downloadHandler(
    filename <- "set_scores.csv",
    
    content = function(file) {
      write.csv(res_ext_dt[, colnames(res_ext_dt) %in% c("RANK (PENETRANCE_NLL)", "RS_IDS", score_names), with=FALSE], 
                file, row.names = FALSE)
    }
  )
  
  # Set scores plot download
  output$set_scores_plt_download <- downloadHandler(
    filename <- "set_scores.png",
    
    content = function(file) {
      png(file,
          width = input$shiny_width,
          height = input$shiny_height)
      
      print(set_scores_plt_re())
      dev.off()
    }
  )
  
  ## Candidate table 
  
  # Copy option for candidates table
  cand_table_prep <- res_ext_dt[, c(1:16)]
  
  cand_table_prep$cp <- vapply(1L:nrow(cand_table_prep), function(i){
    as.character(
      rclipButton(
        paste0("clipbtn_", i), 
        label = "Copy", 
        clipText = cand_table_prep[i, "RS_IDS"], 
        icon = icon("copy", lib = "glyphicon"),
        class = "btn-primary btn-sm"
      )
    )
  }, character(1L))
  
  cand_table_prep <- cbind(cand_table_prep[, c(1:2)], cand_table_prep$cp, cand_table_prep[, c(3:16)])
  colnames(cand_table_prep)[3] <- "."

  
  # Output for results table with candidates
  output$cand_tbl_epi <- DT::renderDT({
    fancyDT(cand_table_prep, "candidate_results")
  })
  
  # Network stats download
  output$cand_tbl_download <- downloadHandler(
    filename <- "candidates.csv",
    
    content = function(file) {
      write.csv(res_ext_dt[, c(1:16)], file, row.names = FALSE)
    }
  )
  


  ### Network visualization ###

  ## Options for networks ##

  
  # rs_ids
  selected_SNPs_net_epi <- reactive({
    
    strsplit(input$selected_SNPs_net_epi, "\\; |\\;|\\, |\\,| ")[[1]]
    
  })

  
  # Priority filter for network centralities 
  centr_filter_names <- reactive({
  
    # check which filters should be used
    centr_filter_merged <- c()
    centr_filter_asc_merged <- c()
    
    if(input$centr_filter_degree != 0){
      centr_filter_merged <- c("degree" = "degree")
      centr_filter_asc_merged <- c("degree" = input$centr_filter_degree)
      }
    
    if(input$centr_filter_eigenvector != 0){
      centr_filter_merged <- c(centr_filter_merged, "eigenvector" = "eigenvector")
      centr_filter_asc_merged <- c(centr_filter_asc_merged, "eigenvector" = input$centr_filter_eigenvector)
    }
    
    if(input$centr_filter_pagerank != 0){
      centr_filter_merged <- c(centr_filter_merged, "pagerank" = "pagerank")
      centr_filter_asc_merged <- c(centr_filter_asc_merged, "pagerank" = input$centr_filter_pagerank)
    }
    
    centr_filter_names <- list(centr_filter_merged, centr_filter_asc_merged)
 
    centr_filter_names
    
  })
  
  # Show only the centralities that are not none 
  output$centr_priority <- renderUI({
    rank_list(input_id = "centr_priority_opt",
              text = "Centrality priority", 
              labels = unlist(centr_filter_names()[[1]]), 
              options = sortable_options(multiDrag = TRUE))
    
  })
  
  # Disable/enable degree_slider 
  observeEvent(input$centr_filter_degree, {
    if(input$centr_filter_degree == 0){
      shinyjs::disable(id = "degree_slider")
    } else {
      shinyjs::enable(id = "degree_slider")
    }
  })
 
  
  # Disable/enable eigenvector_slider 
  observeEvent(input$centr_filter_eigenvector, {
    if(input$centr_filter_eigenvector == 0){
      shinyjs::disable(id = "eigenvector_slider")
    } else {
      shinyjs::enable(id = "eigenvector_slider")
    }
  })
  
  # Disable/enable pagerank_slider 
  observeEvent(input$centr_filter_pagerank, {
    if(input$centr_filter_pagerank == 0){
      shinyjs::disable(id = "pagerank_slider")
    } else {
      shinyjs::enable(id = "pagerank_slider")
    }
  })
  
  
  # Disable/enable number of non-candidate neighbors (num_res)
  centr_filter_check <- reactive({
    list(input$centr_filter_degree,
         input$centr_filter_eigenvector, 
         input$centr_filter_pagerank)
  })
  
  observeEvent(centr_filter_check(), {
    if((input$centr_filter_degree == 0) &
       (input$centr_filter_eigenvector == 0) &
       (input$centr_filter_pagerank == 0)){
      updateSliderInput(session, inputId = "num_res", value = 0)
      shinyjs::disable(id = "num_res")
    } else {
      shinyjs::enable(id = "num_res")
    }
  })

  
  # Network preparation
  network_re_epi <- eventReactive(input$snp_net_action, {

    SNPs_to_validate <- strsplit(input$selected_SNPs_net_epi, "\\; |\\;|\\, |\\,| ")[[1]]
    
    # Validate Input
    if(sum(!grepl("^rs\\d*$", SNPs_to_validate)) > 0){
      updateTextInput(inputId = "selected_SNPs_net_epi", 
                      value = "",
                      placeholder = "rs7596121;rs41440544;rs41501252;rs16976638;rs16976644;rs10518828;rs16976648"
      )
      shinyalert(title = "Incorrect input.", 
                 text = "Please enter rsIDs in the correct input format.", 
                 type = "info", 
                 closeOnEsc = TRUE
      )
  
      # if input is incorrect show no network
      show_network_from_id(nodes, edges, c(""), filter = FALSE,
                           centr_filter = c("degree"),
                           ascending = c(-1), num_res = 0)
      
    } else if(sum(! gsub("rs","", SNPs_to_validate) %in% nodes[,rs_id]) > 0) {
      updateTextInput(inputId = "selected_SNPs_net_epi", 
                      value = "",
                      placeholder = "rs7596121;rs41440544;rs41501252;rs16976638;rs16976644;rs10518828;rs16976648"
      )
      shinyalert(title = "Invalid rsID/rsIDs.", 
                 text = "Please enter valid rsIDs, which exist in the network, as input above.", 
                 type = "info", 
                 closeOnEsc = TRUE
      )
      
      # if input is incorrect show no network
      show_network_from_id(nodes, edges, c(), filter = TRUE,
                           centr_filter = c("degree"),
                           ascending = c(-1), num_res = 0)
      
      
    } else {
  
      # calculate closeness for all SNPs in the current candidate set (if not already calculated)
      
      # remove rs from rsIDs
      clo_cands <- gsub("rs", "", selected_SNPs_net_epi())
      
      # candidates for which the closeness has to be calculated
      clo_cands <- nodes[rs_id %in% clo_cands & is.na(closeness), rs_id]
      
      # calculate closeness for closeness candidates 
      nodes <<- get_cand_centr(i_graph, nodes, clo_cands, centr_measure = "closeness", mode = "all")
      nodes$closeness.x <<- nodes[, ifelse(is.na(closeness.x), closeness.y, closeness.x)]
      nodes$closeness.y <<- NULL
      colnames(nodes)[6] <<- c("closeness")
      
    
    # get into for network
    net_info <- tryCatch({

      if((input$centr_filter_degree == 0) &
         (input$centr_filter_eigenvector == 0) &
         (input$centr_filter_pagerank == 0)){

        # if no filter is selected, show only candidates
        show_network_from_id(nodes, edges, selected_SNPs_net_epi(), filter = TRUE,
                             centr_filter = c("degree"),
                             ascending = c(-1), num_res = 0)
      } else {
        
        nodes_for_priority <- nodes
        

        # select SNPs to keep so they are not filtered out by the centralities
        candidate_snps_to_keep <- selected_SNPs_net_epi()
        candidate_snps_to_keep <- gsub("rs", "", candidate_snps_to_keep)
        nodes_from_cands <- nodes[rs_id %in% candidate_snps_to_keep]
        
        
        if("degree" %in% centr_filter_names()[[1]]){
          # Nodes to prioritize 
          nodes_for_priority <- nodes_for_priority[(degree %between% c(input$degree_slider[1], input$degree_slider[2]))]
        }
        
        if("eigenvector" %in% centr_filter_names()[[1]]){
          # Nodes to prioritize 
          nodes_for_priority <- nodes_for_priority[(eigenvector %between% c(input$eigenvector_slider[1], input$eigenvector_slider[2]))]
        }
        
        if("pagerank" %in% centr_filter_names()[[1]]){
          # Nodes to prioritize 
          nodes_for_priority <- nodes_for_priority[(pagerank %between% c(input$pagerank_slider[1], input$pagerank_slider[2]))]
        }
        
        nodes_for_priority <- unique(rbind(nodes_for_priority, nodes_from_cands))
      
        
        show_network_from_id(nodes_for_priority, edges, selected_SNPs_net_epi(), filter = TRUE,
                             centr_filter = centr_filter_names()[[1]][order(input$centr_priority_opt)],
                             ascending = centr_filter_names()[[2]][order(input$centr_priority_opt)], 
                             num_res = input$num_res
        )

      }


    }, error = function(x) {
      NULL
    })

    net_info
    }
  })


  # Network
  output$network_epi <- renderVisNetwork({
    network_re_epi()$net %>%
      visEvents(click = "function(nodes){
                  Shiny.setInputValue('click_snp', nodes.nodes[0]);
                  ;}"
      )
  })
  
  
  # Download image of snp-snp interaction network for selected candidate SNPs
  output$snp_net_download <- downloadHandler(
    filename <- "snp_net.html",
    
    content = function(con) {
      network_re_epi()$net %>% visSave(con)
      
    }
  )

  # Protein-protein interaction network preparation
  gene_network_re_epi <- reactive({

    net <- tryCatch({

      show_gene_network_from_ids(rs_ids = network_re_epi()$nodes$label,
                                 dbSNP_dt = dbSNP_dt,
                                 biogrid_dt = biogrid_dt)
      
    }, error = function(x) {
      NULL
    })
    
    net
  })
  
  # Download image of snp-snp interaction network for selected candidate SNPs
  output$gene_net_download <- downloadHandler(
    filename <- "gene_net.html",
    
    content = function(con) {
      gene_network_re_epi()$net %>% visSave(con)
      
    }
  )
  

  # Table for current SNP
  current_SNP_tbl_re_epi <- reactive({

    SNP_tbl <- tryCatch({
      # get nodes data table for current SNP
      current_nodes_dt <- nodes[id == input$click_snp]
      
      # remove id column
      current_nodes_dt <- current_nodes_dt[, id := NULL]
      
      # add gene
      current_gene_for_snp <- gene_network_re_epi()[[2]][gene_network_re_epi()[[2]]$rs_id == current_nodes_dt$rs_id, `gene/s`]
      
      # add rs from rs_id
      current_nodes_dt$rs_id <- paste0("rs", current_nodes_dt$rs_id)
      
      # save current rs_id
      current_rs_id <- current_nodes_dt$rs_id
      
      # delete rs_id from this data table for correct order when binding
      current_nodes_dt <- current_nodes_dt[, rs_id := NULL]

      #  table with properties and values
      current_dt <- data.table(property = colnames(current_nodes_dt),
                     value = c(current_nodes_dt))
      
      
      current_dt <- rbind(data.table(property = c("rs_id"), 
                                     value = current_rs_id), 
                          data.table(property = c("gene/s"), 
                                                 value = eval(current_gene_for_snp)), 
                          current_dt)
      
      current_dt[!is.na(value)]
      
    }, error = function(x) {
      data.table()
    })

    SNP_tbl
  })


  # Show table for current SNP
  output$snp_tbl_epi <- renderText({
    as.character(
      div(
        div(
          class="row",
          table_custom(current_SNP_tbl_re_epi(), title = "", 
                       class="col-md-12", show_rownames_val = FALSE),
        )
      )
    )
  })


  # Network stats for current network
  output$snp_net_stats_tbl_epi <- renderText({
    as.character(
      div(
        class="row",
        table_custom(get_net_stats(network_re_epi()$nodes,
                                   data.table(network_re_epi()$edges)[from %in% network_re_epi()$nodes$id & to %in% network_re_epi()$nodes$id]),
                     title = "", class="col-md-12", 
                     show_rownames_val = FALSE),
      )
    )
  })


  # Network stats download
  output$snp_net_stats_download <- downloadHandler(
    filename <- "snp_network_stats.csv",
    
    content = function(file) {
      write.csv(get_net_stats(network_re_epi()$nodes,
                              data.table(network_re_epi()$edges)[from %in% network_re_epi()$nodes$id & to %in% network_re_epi()$nodes$id]), file, row.names = FALSE)
    }
  )

  ### Databases (for current SNP as input)
  ## SNP databases

  # Current rs_id of SNP
  current_SNP_re_epi <- reactive({

    current_snp <- tryCatch({

      current_nodes_dt <- nodes[id == input$click_snp]
      current_nodes_dt$rs_id <- paste0("rs", current_nodes_dt$rs_id)

      current_rs_id <- current_nodes_dt$rs_id
      current_rs_id

    }, error = function(x) {
      ""
    })

    current_snp
  })



  # dbSNP
  output$dbSNP_link_epi <- renderUI({
    tags$a(
      href=ifelse(current_SNP_re_epi() == "",
                  "https://www.ncbi.nlm.nih.gov/snp/",
                  paste0("https://www.ncbi.nlm.nih.gov/snp/?term=",
                         current_SNP_re_epi())),
      target="_blank",
      tags$img(src = dbSNP_img_path,
               alt = "dbSNP",
               style = "width:35%")
    )
  })

  # ClinVar
  output$clinvar_link_epi <- renderUI({
    tags$a(
      href=ifelse(current_SNP_re_epi() == "",
                  "https://www.ncbi.nlm.nih.gov/clinvar/",
                  paste0("https://www.ncbi.nlm.nih.gov/clinvar/?term=",
                         current_SNP_re_epi())),
      target="_blank",
      tags$img(src = clinvar_img_path,
               alt = "ClinVar",
               style = "width:35%")
    )
  })

  
  # SNPmap
  output$snpmap_link_epi <- renderUI({
    tags$a(
      href=ifelse(current_SNP_re_epi() == "",
                  "http://snp.nbscn.org/",
                  paste0("http://snp.nbscn.org/rs/id/",
                         gsub("rs","",current_SNP_re_epi()))),
      target="_blank",
      tags$img(src = snpmap_img_path,
               alt = "SNPmap",
               style = "width:35%"
               )
    )
  })
  

  ### Information on selected SNPs
  
  output$info_on_selected_SNPs_box <- renderUI({
    
    
    SNPs_to_validate <- strsplit(input$selected_SNPs_net_epi, "\\; |\\;|\\, |\\,| ")[[1]]
    
    # Validate Input
    if(input$selected_SNPs_net_epi == "" | sum(!grepl("^rs\\d*$", SNPs_to_validate)) > 0){
      div(
        h5("Please enter rsIDs in the correct input format.")
      )
    } else if(sum(! gsub("rs","", SNPs_to_validate) %in% nodes[,rs_id]) > 0) {
      div(
        h5("Please enter valid rsIDs, which exist in the network, as input above.")
      )
    } else if(input$snp_net_action == 0) {
      div(
        h5("Please click on 'Create network'.")
      )
    } else {
  
  div(
    h3("Overview of selected SNPs"),
    h4(align = "justify",
       style = "font-size: 500; line-height: 1.5;",
       HTML(
         "To obtain an overview of all selected candidate SNPs and their filtered 
         neighborhood, the following table presents the same properties as the 
         one for the single selected SNPs. However, this table includes all 
         SNPs of the current SSI network."
       )
    ),
    
    # Download button for network statistics
    downloadBttn("rs_id_gene_tbl_download", 
                 "Download", 
                 style = "unite",
                 color = "primary",
                 size = "sm",
                 block = FALSE,
                 no_outline = TRUE),
    
    # Table with rs_id gene mapping and scores
    shinycssloaders::withSpinner(
      DTOutput("rs_id_gene_tbl_epi")
    ), 
    
    # Plot with marked SNPs in whole SNP-SNP interaction network 
    # distribution for given centrality 
    # Select box to select the score to display in the snp_scores_plt
    h3("Centrality distribution with selected SNPs"),
    
    
    h4(align = "justify",
       style = "font-size: 500; line-height: 1.5;",
       HTML(
         "In the following the four different centrality distributions of the 
         entire SSI network can be selected. However, in this case the candidate 
         SNPs are highlighted according to the respective centrality to 
         contextualize the centralities of these SNPs."
       )
    ),
    
    fluidRow(
      column(3, 
             selectInput("centr_select", label = h5("Select the centrality to display"), 
                         choices = list("Degree" = "degree",
                                        "Eigenvector" = "eigenvector",
                                        "Pagerank" = "pagerank", 
                                        "Closeness" = "closeness"
                         ), 
                         selected = "degree")
             )
    ),
 
    
    h5(' '), 
    
    # Download button for centrality distribution data with candidate SNPs
    downloadBttn("centr_dist_cand_data_download", 
                 "Download data", 
                 style = "unite",
                 color = "primary",
                 size = "sm",
                 block = FALSE,
                 no_outline = TRUE),
    
    h5(' '), 
    
    # Download button for centrality distribution plots with candidate SNPs
    downloadBttn("centr_dist_cand_plt_download", 
                 "Download plot", 
                 style = "unite",
                 color = "primary",
                 size = "sm",
                 block = FALSE,
                 no_outline = TRUE),
    
    
    h5(' '), 
    
    
    shinycssloaders::withSpinner(
      plotOutput(outputId = "cand_centr_plt", height = "500px")
    ),
    
    
    # Select box to select the score to display in the snp_scores_plt
    h3("Statistical scores"),
    
    h4(align = "justify",
       style = "font-size: 500; line-height: 1.5;",
       HTML(
         "After selecting one of the 14 statistical scores, the score 
         distribution of the selected candidate SNPs will be displayed with 
         respect to the considered score. The x-axis of each plot contains the 
         rsIDs of the candidate SNPs, and the selected score is shown on the y-axis. 
         By hovering over a data point, additional information can be retrieved,
         including the rsID, the gene or genes impacted by the SNP, and the
         current selected score for the SNP."
       )
    ),
    
    selectInput("score_select", label = h5("Select the score to display"), 
                choices = list("PENETRANCE_NLL" = "PENETRANCE_NLL",
                               "PENETRANCE_LLH" = "PENETRANCE_LLH",
                               "PENETRANCE_AIC" = "PENETRANCE_AIC",
                               "PENETRANCE_BIC" = "PENETRANCE_BIC",    
                               "REGRESSION_NLL" = "REGRESSION_NLL",
                               "REGRESSION_LLH" = "REGRESSION_LLH", 
                               "REGRESSION_AIC" = "REGRESSION_AIC",
                               "REGRESSION_BIC" = "REGRESSION_BIC",
                               "REGRESSION_NLL_GAIN" = "REGRESSION_NLL_GAIN",
                               "REGRESSION_LLH_GAIN" = "REGRESSION_LLH_GAIN",
                               "REGRESSION_AIC_GAIN" = "REGRESSION_AIC_GAIN",
                               "REGRESSION_BIC_GAIN" = "REGRESSION_BIC_GAIN",
                               "BAYESIAN" = "BAYESIAN",
                               "VARIANCE (CHI-SQUARE)" = "VARIANCE"), 
                selected = "PENETRANCE_NLL"),
    
    h5(' '), 
    
    # Download button for centrality distribution data with candidate SNPs
    downloadBttn("snp_scores_cand_data_download", 
                 "Download data", 
                 style = "unite",
                 color = "primary",
                 size = "sm",
                 block = FALSE,
                 no_outline = TRUE),
    
    h5(' '), 
    
    # Download button for centrality distribution plots with candidate SNPs
    downloadBttn("snp_scores_cand_plt_download", 
                 "Download plot", 
                 style = "unite",
                 color = "primary",
                 size = "sm",
                 block = FALSE,
                 no_outline = TRUE),
    
    h5(' '), 
    
    shinycssloaders::withSpinner(
      plotlyOutput("snp_scores_plt")
    ),
    
    
  )
    
  }
    
  })
  
  
 
  # rs id gene table
  rs_id_gene_table <- reactive({
    rs_id_gene_table <- gene_network_re_epi()[[2]]
    rs_id_gene_table <- merge(rs_id_gene_table, nodes, by="rs_id", all.x = TRUE)
    rs_id_gene_table$rs_id <- paste0("rs", rs_id_gene_table$rs_id)
    rs_id_gene_table <- rs_id_gene_table[, id := NULL]
  })
  
  
  
  # Table with rs_id and corresponding gene/s, centralities and scores
  # for the current SNP-SNP interaction network
  output$rs_id_gene_tbl_epi <- DT::renderDataTable({
    fancyDT(rs_id_gene_table(), "rs_id_gene_centr_scores")
  })

  # Network stats download
  output$rs_id_gene_tbl_download <- downloadHandler(
    filename <- "selected_candidates_info.csv",
    
    content = function(file) {
      write.csv(rs_id_gene_table(), file, row.names = FALSE)
    }
  )
  
  
  # Plot with marked SNPs in whole SNP-SNP interaction network 
  # distribution for given centrality 
  cand_centr_plt_re <- reactive({
    plot_centr_with_cand(nodes = nodes,
                         rs_ids = selected_SNPs_net_epi(),
                         centrality = input$centr_select,
                         num_bins = 50)
  })
  
  # Distribution of selected centrality in the SSI network with candidate SNPs
  output$cand_centr_plt <- renderPlot({
    cand_centr_plt_re()
  })
  
  
  # Download for distribution of selected centrality data with candidate SNPs
  output$centr_dist_cand_data_download <- downloadHandler(
    filename = function() { paste(input$centr_select, '_centr_dist_cand.csv', sep='') },
    
    content = function(file) {
      tmp_nodes <- nodes[, .(rs_id, get(input$centr_select))]
      colnames(tmp_nodes) <- c('rsID', input$centr_select)
      write.csv(tmp_nodes, file, row.names = FALSE)
    }
  )
  
  # Download for distribution of selected centrality plot with candidate SNPs
  output$centr_dist_cand_plt_download <- downloadHandler(
    filename = function() { paste(input$centr_select, '_centr_dist_cand.png', sep='') },
    
    content = function(file) {
      png(file,
          width = input$shiny_width,
          height = input$shiny_height)
      
      print(cand_centr_plt_re())
      dev.off()
    }
  )
  
  ## Plot with statistical scores (selection) of the individual candidate SNPs
  output$snp_scores_plt <- renderPlotly({
    rs_id_score_gene_dt <- rs_id_gene_table()
    
    rs_id_score_gene_dt$rs_id <- gsub("rs", "", rs_id_score_gene_dt$rs_id)
    
    plot_snp_scores(rs_ids = selected_SNPs_net_epi(), 
                    snp_scores_dt = rs_id_score_gene_dt, 
                    score_name = input$score_select, 
                    show_gene = TRUE) 
  })
  
  
  # Download for statistical score data of selected scores with candidate SNPs
  output$snp_scores_cand_data_download <- downloadHandler(
    filename = function() { paste(input$score_select, '_score_cands.csv', sep='') },
    
    content = function(file) {
      tmp_tbl <- rs_id_gene_table()[, .(rs_id, `gene/s` , get(input$score_select))]
      colnames(tmp_tbl) <- c('rsID', 'gene/s', input$score_select)  
      write.csv(tmp_tbl, file, row.names = FALSE)
    }
  )
  
  # Statistical score plot static for download
  snp_scores_plt_re <- reactive({
    rs_id_score_gene_dt <- rs_id_gene_table()
    
    rs_id_score_gene_dt$rs_id <- gsub("rs", "", rs_id_score_gene_dt$rs_id)
    
    plot_snp_scores_static(rs_ids = selected_SNPs_net_epi(),
                           snp_scores_dt = rs_id_score_gene_dt, 
                           score_name = input$score_select)
  })
  

  # Download for statistical score plot of selected scores with candidate SNPs
  output$snp_scores_cand_plt_download <- downloadHandler(
    filename = function() { paste(input$score_select, '_score_cands.png', sep='') },
    
    content = function(file) {
      png(file,
          width = input$shiny_width,
          height = input$shiny_height)
     
      print(snp_scores_plt_re())
      
      dev.off()
    }
  )
  
  
  ####################################################
  # EDA: PART 3: PROTEIN-PROTEIN INTERACTION NETWORK # 
  ####################################################
  
  ## Protein-protein interaction network
  output$gene_gene_box <- renderUI({
    SNPs_to_validate <- strsplit(input$selected_SNPs_net_epi, "\\; |\\;|\\, |\\,| ")[[1]]
    
    # Validate Input
    if(input$selected_SNPs_net_epi == "" | sum(!grepl("^rs\\d*$", SNPs_to_validate)) > 0){
      div(
        h5("Please enter rsIDs in the correct input format.")
      )
    } else if(sum(! gsub("rs","", SNPs_to_validate) %in% nodes[,rs_id]) > 0) {
      div(
        h5("Please enter valid rsIDs, which exist in the network, as input above.")
      )
    } else if(input$snp_net_action == 0) {
      div(
        h5("Please click on 'Create network'.")
      )
    } else {
      
      div(
        
        h4(align = "justify",
           style = "font-size: 500; line-height: 1.5;",
           HTML(
             'To construct the PPI network all SNPs of the current SSI network
             were first mapped to genes via dbSNP, and the edges were introduced 
             using the PPI network from BioGRID. In the PPI network, the nodes 
             represent the proteins, and they are labelled accordingly. <br>
             By using the "Select by id" field in the top left corner of the 
             network, it is possible to choose a gene of interest, which will
             then be highlighted in the network.'
           )
        ),
        
        shinycssloaders::withSpinner(
          visNetworkOutput(outputId = "gene_network_epi")
        ),
        
        h5(' '), 
        
        # Download button for centrality distribution data with candidate SNPs
        downloadBttn("gene_net_download", 
                    "Download network", 
                    style = "unite",
                    color = "primary",
                    size = "sm",
                    block = FALSE,
                    no_outline = TRUE),
        
        h5(' '),
        
        fluidRow(
          column(4, 
                 # Gene databases
                 box(
                   width = 12, 
                   title = "Gene databases",
                   solidHeader = TRUE, 
                   # background = "red"
                   status = "primary",
                   collapsible = TRUE,
                   collapsed = FALSE,
                   
                   h4(align = "justify",
                      style = "font-size: 500; line-height: 1.5;",
                      HTML(
                        "The following four databases allow for the 
                        analysis of the currently viewed gene. 
                        If no node is selected, the home page 
                        will be displayed."
                      )
                   ),
                   
                   fluidRow(
                     column(6, 
                            # Genecards
                            uiOutput(outputId = "genecards_link_epi")
                     ), 
                     h4(""),
                     column(6, 
                            # Cytoscape
                            uiOutput(outputId = "cytoscape_link_epi")
                     )
                   ), 
                   fluidRow(
                     column(6, 
                            # Genbank
                            uiOutput(outputId = "genbank_link_epi")
                     ), 
                     h4(""),
                     column(6, 
                            # Gene ontology
                            uiOutput(outputId = "geneontology_link_epi")
                     )
                   )
                 )
          ), 
          column(4,
                 # Gene enrichment
                 box(
                   width = 12, 
                   title = "Gene enrichment",
                   solidHeader = TRUE, 
                   # background = "red"
                   status = "primary",
                   collapsible = TRUE,
                   collapsed = FALSE,
                   
                   h4(align = "justify",
                      style = "font-size: 500; line-height: 1.5;",
                      HTML(
                        "The following database allows for performing gene set 
                        enrichment analysis on the currently affected genes."
                      )
                   ),
                   
                   # gprofiler
                   uiOutput(outputId = "gprofiler_link_epi")
                   
                 )
          ), 
          column(4, 
                 # SNP information e.g., SNP-gene mapping
                 box(
                   width = 12, 
                   title = "SNP information",
                   solidHeader = TRUE, 
                   # background = "red"
                   status = "primary",
                   collapsible = TRUE,
                   collapsed = FALSE,
                   
                   h4(align = "justify",
                      style = "font-size: 500; line-height: 1.5;",
                      HTML(
                        "Enter an rsID to retrieve the corresponding gene or genes."
                      )
                   ),
                   
                   textInput(inputId = "selected_SNP_bio", 
                             label = "Enter SNP:", 
                             width = "100%",
                             placeholder = "rs35664415"),
                   
                   
                   div("Corresponding gene(s):", 
                       shinycssloaders::withSpinner(
                         htmlOutput("snp_info_epi")
                       ))
                   
                   )
                 )
        )
      )
    }
    
  })
  
  
  # Protein-protein interaction network with genes from candidate SNPs
  output$gene_network_epi <- renderVisNetwork({
    gene_network_re_epi()[[1]] %>%
      visEvents(click = "function(nodes){
                  Shiny.setInputValue('click_gene', nodes.nodes[0]);
                  ;}"
      )
  })

  # Table with gene/protein and id of the protein-protein interaction network
  gene_network_nodes_re_epi <- reactive({

    gene_id_nodes <- tryCatch({
      gene_nodes <- gene_network_re_epi()[[2]]$`gene/s`
      gene_nodes <- data.table(gene = str_sort(unique(gene_nodes)))
      gene_nodes[, id := seq(0, dim(gene_nodes)[1]-1)]
      gene_nodes

    }, error = function(x) {
      NULL
    })

    gene_id_nodes
  })
  

  
  # Current gene name
  current_gene_re_epi <- reactive({

    current_gene <- tryCatch({
      gene_network_nodes_re_epi()[id == input$click_gene, gene]

    }, error = function(x) {
      ""
    })

    current_gene
  })
  
  
  
  

  ## Gene databases

  # GeneCards
  output$genecards_link_epi <- renderUI({
    tags$a(
      href=ifelse(current_gene_re_epi() == "",
                  "https://www.genecards.org/",
                  paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                         current_gene_re_epi())),
      target="_blank",
      tags$img(src = genecards_img_path,
               alt = "GeneCards",
               style = "width:100%")
    )
  })

  # GenBank
  output$genbank_link_epi <- renderUI({
    tags$a(
      href=ifelse(current_gene_re_epi() == "",
                  "https://www.ncbi.nlm.nih.gov/genbank/",
                  paste0("https://www.ncbi.nlm.nih.gov/gene/?term=",
                         current_gene_re_epi())),
      target="_blank",
      tags$img(src = genbank_img_path,
               alt = "GenBank",
               style = "width:100%")
    )
  })


  # WikiPathways


  # Gene ontology
  output$geneontology_link_epi <- renderUI({
    tags$a(
      href=ifelse(current_gene_re_epi() == "",
                  "http://geneontology.org/",
                  paste0("http://amigo.geneontology.org/amigo/medial_search?q=",
                         current_gene_re_epi(), "&searchtype=all")),
      target="_blank",
      tags$img(src = geneontology_img_path,
               alt = "GeneOntology",
               style = "width:100%")
    )
  })


  # Cytoscape
  output$cytoscape_link_epi <- renderUI({
    tags$a(
      href=ifelse(current_gene_re_epi() == "",
                  "https://cytoscape.org/",
                  paste0("https://www.ndexbio.org/iquery/?genes=",
                         current_gene_re_epi(), "%2C")),
      target="_blank",
      tags$img(src = cytoscape_img_path,
               alt = "Cytoscape",
               style = "width:100%")
    )
  })

  # gprofiler
  output$gprofiler_link_epi <- renderUI({
    tags$a(
      href=ifelse(current_gene_re_epi() == "",
                  "https://biit.cs.ut.ee/gprofiler/page/r",
                  paste0("https://biit.cs.ut.ee/gprofiler/gost?organism=hsapiens&query=", paste(gene_network_nodes_re_epi()$gene, collapse = "%20%0A"), "%0A%0A&ordered=false&all_results=false&no_iea=false&combined=false&measure_underrepresentation=false&domain_scope=annotated&significance_threshold_method=g_SCS&user_threshold=0.05&numeric_namespace=ENTREZGENE_ACC&sources=GO:MF,GO:CC,GO:BP,KEGG,TF,REAC,MIRNA,HPA,CORUM,HP,WP&background=&no_evidences=false")),
      target="_blank",
      tags$img(src = gprofiler_img_path,
               alt = "gprofiler",
               style = "width:70%")
    )
  })
  
  
  # SNP information e.g., SNP-gene mapping
  output$snp_info_epi <- renderUI({
    HTML(
      paste(c(get_gene_from_snp(input$selected_SNP_bio, dbSNP_dt)), "", sep = "<br/>")
    )
  })
  
  
  output$genome_browser_frame <- renderUI({
    tags$iframe(src="https://vps-zap841713-2.zap-srv.com/genome-browser?genome=hg38", 
                height=input$shiny_height, 
                width=input$shiny_width)
  })
  
  


  # Protection so app closes correctly without crashing everything
  session$onSessionEnded(function() {
    stopApp()
  })

}


shinyApp(ui, server)
