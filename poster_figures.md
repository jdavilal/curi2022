## Method Figures

plot1

    #Load the P3, Sig10,6,5 data
    p3 <- read_csv("Data/P3_sig.csv")
    mutational_sig <- read_csv("mutational_sig_oec/supplied_data/mutational.sig.csv")

    p3_withcontext <- mutational_sig %>%
      filter(sample == "s1133903279") %>%
      select(1,4,7) %>%
      mutate(p3 = prob + 0.0001) %>%
      select(1,3,4)
    p3_withcontext <- tibble::rowid_to_column(p3_withcontext, "ID")

    p3 <- tibble::rowid_to_column(p3, "ID")
    p3_final <- full_join(p3, p3_withcontext, by = "ID")
    p3_final <- p3_final %>%
      mutate(p3 = p3.x) %>%
      select(2:4, 7, 9)

    #plot1
    p1 <- ggplot(p3_final, aes(x = sig10, y = p3)) +
      geom_point(alpha = 0.3) + 
      geom_text(data = subset(p3_final, p3 > 0.12 & sig10 < 0.1), aes(label=ext.context), vjust = 2, hjust = -0.25) +
      geom_text(data = subset(p3_final, p3 > 0.14 & sig10 > 0.2), aes(label=ext.context), vjust = 2) +
      geom_text(data = subset(p3_final, sig10 > 0.3), aes(label=ext.context), vjust = 2, hjust = 1) +
      geom_smooth(method = "lm", color = "mediumpurple3", level = 0.5) +
      labs(x = "SBS10",
           y = "Sample Profile") +
      theme_bw() +
      theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
            axis.title = element_text(size = 13)) +
      theme(panel.grid.minor.x = element_line(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_blank())
    p1

![](poster_figures_files/figure-markdown_strict/Method%20Figures-1.png)

    #ggsave("method1.png", device='png', dpi=700)

plot2

    #plot2
    p2 <- ggplot(p3_final, aes(x = sig6, y = p3)) +
      geom_point(alpha = 0.3) +
      geom_text(data = subset(p3_final, p3 > 0.12 & sig6 < 0.1), aes(label=ext.context), vjust = 2) +
      geom_text(data = subset(p3_final, p3 > 0.12 & sig6 > 0.1), aes(label=ext.context), vjust = 2) +
      geom_smooth(method = "lm", color = "mediumpurple3", level = 0.5) +
      labs(x = "SBS6",
           y = "Sample Profile") +
      xlim(0, 0.3) +
      theme_bw() +
      theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
            axis.title = element_text(size = 13)) +
      theme(panel.grid.minor.x = element_line(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_blank())
    p2

![](poster_figures_files/figure-markdown_strict/unnamed-chunk-1-1.png)

    #ggsave("method2.png", device='png', dpi=700)

Method2 figure (plot1 + plot2)

    ggarrange(p1, p2,
                        labels = c("B", "C"),
                        ncol = 2)

![](poster_figures_files/figure-markdown_strict/optional-1.png)

    #ggsave("method1v2.png", device='png', dpi=700)

## Result Figures

Result1 figure

    eval_df <- read_csv("Cancer_Evaluation/Cancer_Evaluation_Data/colorect_cancer_evaluation.csv")

    #Result Figure1
    ggplot(data = eval_df, aes(x = expected_list, y = actual_list, label = sig_list))+
      geom_point(alpha=0.3) +
      geom_smooth(method = "lm", color = "mediumpurple3") +
      labs(x = "Results from Dissolvo", 
            y = "True Proportions",
            title = "Comparison of Dissolvo Output and 
    Sample Colorectal Cancer Profiles") +
      theme_bw() +
      theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
      theme(panel.grid.minor.x = element_line(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_blank())

![](poster_figures_files/figure-markdown_strict/Result%20Figures-1.png)

    #ggsave("result1.png", device='png', dpi=700)

Result2 figure

    sig_levels <- c("SBS1", "SBS5", "SBS15", "SBS40", "SBS44")
    eval_df2 <- eval_df %>%
      mutate(sig_list = factor(sig_list, sig_levels)) %>%
      group_by(sig_list) %>%
      mutate(n = n()) %>%
      filter(n > 65 & !is.na(sig_list))

    #Result Figure2
    ggplot(data = eval_df2, aes(x = expected_list, y = actual_list))+
       geom_point(alpha=0.3) +
       geom_smooth(method = "lm", color = "mediumpurple3") +
       facet_wrap(~sig_list) +
       labs(x = "Results from Dissolvo", y = "True Proportions", title = "Dissolvo Output and Sample Profiles
    of Selected Mutational Signatures") +
      theme_bw() +
      theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
      theme(panel.grid.minor.x = element_line(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_blank())

![](poster_figures_files/figure-markdown_strict/unnamed-chunk-2-1.png)

    #ggsave("result2.png", device='png', dpi=700)

Result3 figure

    acc <- read_csv("Data/top10_cancer_correlation.csv")
    acc <- acc %>%
      mutate(all_cancer = fct_recode(all_cancer,
                                     "Lung" = "Lung_cancer",
                                     "Head" = "Head_cancer",
                                     "Bladder" = "Bladder_cancer", 
                                     "Skin" = "Skin_cancer",
                                     "Liver" = "Liver_cancer",
                                     "Uterus" = "Uterus_cancer",
                                     "Breast" = "Breast_cancer",
                                     "Colorectal" = "Colorectal_cancer", 
                                     "Stomach" = "Stomach_cancer",
                                     "Prostate" = "Prostate_cancer"))


    #Result Figure3
    ggplot(acc, aes(y= all_corr, x= fct_reorder(all_cancer, all_corr))) + 
      geom_bar(stat="identity", color = "black", fill = "#D1C4E9") +
      coord_flip(ylim=c(0.92,1), expand = c(0, 0)) +
      labs(x = "Cancer Type", y = "Correlation Values", title = "Correlations between Dissolvo and 
    True Proportion for Ten Cancer Types") +
      theme_bw() +
      theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
      theme(panel.grid.minor.x = element_line(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_blank())

![](poster_figures_files/figure-markdown_strict/unnamed-chunk-3-1.png)

    #ggsave("result3.png", device='png', dpi=700)
