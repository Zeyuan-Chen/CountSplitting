library("hash")
library("grid")
library("rlang")
library("scales")
library("ggpubr")
library("ggplot2")
library("extrafont")
library("hrbrthemes")

# p1 is a matrix of num.genes.1 by 1, p2 is a matrix of num.genes.2 by 1
# extract shared genes' pvals
pvals_intersect = function(p1, p2){
    p1 = as.matrix(p1)
    p1 = p1[!is.na(p1), ,drop = F]
    
    p2 = as.matrix(p2)
    p2 = p2[!is.na(p2), ,drop = F]
    keep.genes = intersect(rownames(p1), rownames(p2))
    
    return(cbind(p1 = p1[keep.genes, ],
                 p2 = p2[keep.genes, ]))
}

#pvals is a list of pvals
plot_qq <- function(pvals, labels, ggarrange.nrow, ggarrange.ncol, alpha = 0.05){
    qqplots <- lapply(1:length(pvals), function(p){
        significance_th <- alpha/length(pvals[[p]])
        
        df <- data.frame(pvals.obs = -log10(sort(pvals[[p]])), 
                         pvals.exp = -log10(sort((1:length(pvals[[p]]))/length(pvals[[p]]))));
        qqplot <- ggplot(df, aes(x = pvals.exp, y = pvals.obs)) +
                  stat_binhex(geom = "point", bins=1000, size=1) +
                  geom_hline(yintercept=-log10(significance_th), linetype="dashed", color = "red", size=1) + 
                  geom_abline() +  
                  theme_bw() +
                  ggtitle(labels[p]) +
                  theme(plot.title = element_text(hjust = 0.5)) + 
                  guides(fill="none") +
                  xlab(expression(Expected~-log[10](P))) + ylab(expression(Observed~-log[10](P)))
        return(qqplot)
    })
    return(ggarrange(plotlist = qqplots, ncol = ggarrange.ncol, nrow = ggarrange.nrow))
}

#pvals is a list of matrices each matrix has 2 columns, the 2 pvals you want to compare, first col is x, second is y
plot_qq_compare <- function(pvals, labels, ggarrange.nrow, ggarrange.ncol, alpha = 0.05){
    qqplots <- lapply(1:length(pvals), function(p){
        significance_th <- alpha/nrow(pvals[[p]])
        df <- data.frame(pvals.1 = -log10(pvals[[p]][,1]), 
                         pvals.2 = -log10(pvals[[p]][,2]));
        qqplot <- ggplot(df, aes(x = pvals.1, y = pvals.2)) +
                  stat_binhex(geom = "point", bins=1000, size=1) +
                  geom_abline() +
                  geom_hline(yintercept=-log10(significance_th[[1]]), linetype="dashed", color = "red", size=1)+
                  geom_vline(xintercept=-log10(significance_th[[1]]), linetype="dashed", color = "red", size=1)+
                  theme_bw() +
                  ggtitle(labels[p]) +
                  theme(plot.title = element_text(hjust = 0.5)) +
                  guides(fill="none") +
                  xlab(expression(pval1~-log[10](P))) + ylab(expression(pval2~-log[10](P))) + 
                  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 2))+
                  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 2))
        return(qqplot)
    })
    return(ggarrange(plotlist = qqplots, ncol = ggarrange.ncol, nrow = ggarrange.nrow))
}