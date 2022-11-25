source("code/genus_processed.R")
sig_genera <- composite %>%
    nest(data = -taxonomy) %>%
    mutate(test = map(.x=data, ~wilcox.test(rel_abund~srn, dat=.x) %>%
                          tidy)) %>%
    unnest(test) %>%
    mutate(p.adjust = p.adjust(p.value, method = "BH")) %>%
    filter(p.adjust < 0.05) %>%
    select(taxonomy, p.adjust)


composite %>%
    inner_join(sig_genera, by="taxonomy") %>%
    mutate(rel_abund = 100 * (rel_abund + 1/20000),
           taxonomy = str_replace(taxonomy, "(.*)", "*\\1*"),
           taxonomy = str_replace(taxonomy, "\\*(.*)_unclassified\\*",
                                  "Unclassified<br>*\\1*"),
           srn = factor(srn, levels = c(T, F))) %>%
    ggplot(aes(x = rel_abund, y =taxonomy, color = srn, fill=srn)) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.8,
                                                jitter.width = 0.3),
                shape=21) +
    stat_summary(fun.data = median_hilow, fun.args = list(conf.int=0.5),
                 geom = "pointrange",
                 position = position_dodge(width=0.8),
                 color= "black",
                 show.legend = FALSE) +
    scale_x_log10() +
    scale_color_manual(NULL,
                       breaks = c(F, T),
                       values = c("gray", "dodgerblue"),
                       labels = c("Healthy", "SRN")) +
    scale_fill_manual(NULL,
                       breaks = c(F, T),
                       values = c("gray", "dodgerblue"),
                       labels = c("Healthy", "SRN")) +
    labs(x = "Relative Abundance (%)", y=NULL) +
    theme_classic() +
    theme(
        axis.text.y = element_markdown()
    )


ggsave("figures/significant_genera.tiff", width = 6, height = 4)
    
 



get_sens_spec <- function(threshold, score, actual, direction) {
    
    # threshold <- 100
    # score <- test$score
    # actual <- test$srn
    # direction <- "greaterthan"
    
    predicted <- if(direction == "greaterthan") {
        score > threshold
    } else{
        score < threshold
    }
    
    tp <- sum(predicted & actual)
    tn <- sum(!predicted & !actual)
    fp <- sum(predicted & !actual)
    fn <- sum(!predicted & actual)
    
    specifity <- tn / (tn + fp)
    sensitivity <- tp / (tp + fn)
    
    tibble("specifity" = specifity,
           "sensitivity" = sensitivity)
    
}

    
roc_data <- composite %>%
    inner_join(sig_genera, by="taxonomy") %>%
    select(group, taxonomy, rel_abund, fit_result, srn) %>%
    pivot_wider(names_from = taxonomy, values_from = rel_abund) %>%
    pivot_longer(cols = -c(group, srn), names_to = "metric",
                 values_to = "score") %>%
    # filter(metric == "fit_result") %>%
    nest(data = -metric) %>%
    mutate(direction = if_else(metric == "Lachnospiraceae_unclassified",
                               "lessthan", "greaterthan")) %>%
    mutate(roc_data = map2(.x=data, .y=direction, 
                           ~get_roc_data(.x, .y))) %>%
    unnest(roc_data) %>%
    select(metric, specifity, sensitivity)

    
get_roc_data <- function(x, direction) {
    
    # x <- test
    # direction <-  "greaterthan"
    
    
    threshold <- unique(x$score) %>% 
        sort()
    
    map_dfr(.x=threshold, ~get_sens_spec(.x, x$score, x$srn, direction)) %>%
        rbind(c(specifity=0, sensitivity=1))
    
}


# get_sens_spec(100, test$score, test$srn, "greaterthan")
# get_roc_data(test, "greaterthan")


roc_data %>%
    ggplot(aes(x=1-specifity, y=sensitivity, color=metric)) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, color="grey") +
    theme_classic()

ggsave("figures/roc_figures.tiff", width = 6, height = 4)






