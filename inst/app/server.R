#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(magrittr)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

    load(url("https://github.com/brucemoran/bd2han_validation/raw/master/data/BD2HAN_2.DE_ready.RData"))
    
    #load("../../data/BD2HAN_2_RNAseqR/RData/BD2HAN_2.DE_ready.RData")

    agg_log2tpm_ninf_tb <- agg_log2tpm_tb[!is.infinite(rowSums(agg_log2tpm_tb[,3:dim(agg_log2tpm_tb)[2]])),]

    dataset <- shiny::reactive(log2tpm_ninf_tb)

    observeEvent(input$find_gene, ignoreNULL = FALSE, {
        print(paste0("Searching for: ", input$gene))
        if(length(agg_log2tpm_ninf_tb$external_gene_name[input$gene])>0){
            print(paste0("gene_search is: ", input$gene))
        } else {
            shinyalert::shinyalert(paste0("Gene: ", input$gene, " does not exist in the data, please try again"),
                                   type = "warning",
                                   showConfirmButton = TRUE)
        }
    })

    output$genetable <- DT::renderDT({

        agg_log2tpm_ninf_tb

    })
    output$downloadData <- downloadHandler(
      filename = function() {
        paste0('bd2han_validation.', Sys.Date(), '.csv')
      },
      content = function(con) {
        write.csv(agg_log2tpm_ninf_tb, con)
    })
    output$metatable <- DT::renderDT({

        metadata_tb

    })
    output$plotOut <- renderPlot({
        igene <- input$gene

        gene_tb <- dplyr::filter(.data = agg_log2tpm_ninf_tb,
                                 external_gene_name %in% !!igene)
        tgene_df <- data.frame(value = t(gene_tb[,-c(1,2)]),
                           sample = rownames(t(gene_tb[,-c(1,2)]))) %>% dplyr::left_join(., metadata_tb)
        levs <- c("N_STD", "N_LG", "T_STD", "T_LG")
        tgene_df$group_drug <- factor(tgene_df$TISSUE_ARM,
                                  levels = levs)
        ggplot2::ggplot(data = tgene_df, ggplot2::aes(x = TISSUE_ARM, y = value)) +
        ggplot2::geom_boxplot( ggplot2::aes(colour = TISSUE_ARM)) +
        ggsignif::geom_signif(comparisons = list(c("N_STD","N_LG"),
                                                 c("N_STD","T_STD"),
                                                 c("N_STD","T_LG"),
                                                 c("N_LG","T_STD"),
                                                 c("N_LG","T_LG"),
                                                 c("T_STD","T_LG")), map_signif_level=TRUE, method="wilcox", step_increase=0.1) +
        ggplot2::geom_jitter(position =  ggplot2::position_dodge(0.8)) +
        ggplot2::labs(y = "log2TPM", title = paste0("BD2HAN Validation: Expr. per TISSUE_ARM in ", igene)) +
        ggplot2::theme(axis.text.x =  ggplot2::element_text(angle = 45, hjust = 1))
    })
})
