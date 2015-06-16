#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2014-2015 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "mi.gruenstaeudl@gmail.com"
#version = "2015.06.16.1700"


CommT.plotcolors = function (n=6, h=c(0,360)+15) {
    if ((diff(h)%%360) < 1) {
        h[2] = h[2]-360/n
    }
    hcl(h = (seq(h[1], h[2], length=n)), c=100, l=65)
}


CommT.kfdist = function (post_gt_distrs_BEAST,
                         post_gt_distrs_starBEAST,
                         outlier_num = 1,
                         treedist_select = 2) {
  # Descr:  calculate the tree distances
  # Deps:   phangorn::treedist
  #         reshape::melt
  # I/p:    post_gt_distrs_BEAST = object of class phylo; multiple posterior distributions of gene trees from BEAST
  #         post_gt_distrs_starBEAST = object of class phylo; multiple posterior distributions of gene trees from starBEAST
  #         outlier_num
  #         treedist_select = number of the tree distance selected (2 = Kuhner-Felsenstein distance)
  # O/p:    matrix of floats

  n_genes = length(post_gt_distrs_BEAST)
  # 1. Check if valid input
    if (n_genes > 1) {
  # 2. Initialize outmatrix
        out_m = matrix(nrow=length(post_gt_distrs_BEAST[[1]]),
                       ncol=n_genes,
                       NA)
  # 3. Calculate tree distance for each tree pair
        for (h in 1:n_genes) {
            cat("\n\nANALYZING LOCUS '", names(post_gt_distrs_BEAST)[h], "'\n", sep="")
            cat("\nComparing tree pair:\n ")
            for (i in 1:length(post_gt_distrs_BEAST[[h]])) {
                cat(i, " ", sep="")
                out_m[i,h] = phangorn::treedist(post_gt_distrs_BEAST[[h]][[i]],
                                       post_gt_distrs_starBEAST[[h]][[i]],
                                       check.labels = TRUE)[treedist_select]
            }
        }

  # 4. Rename outmatrix
        clmnames = lapply(sprintf("%03d", 1:n_genes), function(x){paste("gene", x, sep="")})
        colnames(out_m) = clmnames

    } else {stop(cat("\n\nERROR: List of posterior distributions of gene trees from BEAST:  not greater than 1.\n"))}

  # 5. Melt outmatrix
    out_df = reshape::melt(out_m, id.vars=c(clmnames))

  # 6. Add grouping variable
    out_df[,ncol(out_df)+1] = rep("regular", length(out_df[,3]))
    outlier_id = paste("gene", sprintf("%03d", as.integer(outlier_num)), sep="")
    out_df[which(out_df[,2]==outlier_id),4] = "outlier"
    colnames(out_df) = c('post_gen', 'gene_id', 'KF_dist', 'grouping_var')
    
  # 7. Return out_df
    return(out_df)
}


CommT.anova = function (in_df) {
  # Descr:  calculate the anova
  # I/p:    in_df = data frame
  # O/p:    a table

    out_list = list()
    for (i in sprintf("%02d", 1:length(unique(in_df[,'gene_id'])))) {
        handle = in_df
        target = paste("gene0", i, sep="")
        handle[which(handle[,2]==target),4] = target
        handle[which(handle[,2]!=target),4] = "blocked"
        aov_results = summary(aov(KF_dist ~ grouping_var + Error(gene_id), data = handle))
        out_list[[target]] = aov_results$'Error: gene'[[1]][1,5]
    }

    out_m = t(as.data.frame(out_list))
    out_m = round(out_m, digits=4)
    colnames(out_m) = "Pr(>F)"

    return(out_m)
}


CommT.legendpos = function (in_data) {
  # Descr:  generate location coordinates for legend
  # I/p:    in_data
  # O/p:    a list
  
    annot_x_list = annot_y_list = xlim_thres_list = c()
    for (i in split(in_data, in_data$gene_id)) {
        xlim_thres_list = c(annot_x_list, quantile(i[,"KF_dist"], probs = c(0.9999)))
        annot_y_list = c(annot_y_list, max(table(cut(i[,"KF_dist"], breaks=100))))
    }
    xlim_thres_pos = quantile(xlim_thres_list, probs = c(0.95))
    annot_x_pos = quantile(xlim_thres_list, probs = c(0.80))[[1]]
    annot_y_pos = max(annot_y_list)
    
    out_l = list("xlim_thres_pos"=xlim_thres_pos, "annot_x_pos"=annot_x_pos, "annot_y_pos"=annot_y_pos)

    return(out_l)
}


CommT.viz = function (in_df, title_str="a_project_name_here", alpha=0.05, legend_text, legend_pos) {
  # Descr:  visualize the tree distances
  # Deps:   ggplot2::ggplot
  # I/p:    in_df = input dataframe
  #         title_str = string with title name
  #         alpha = significance level
  #         legend_text = input legend
  #         legend_pos = legend position
  # O/p:    a plot

  # 0. Parse legend position information
    annot_x_pos = legend_pos$annot_x_pos
    annot_y_pos = legend_pos$annot_y_pos
    xlim_thres = legend_pos$xlim_thres_pos

  # 1. Define colors
    color_specs = CommT.plotcolors(n=2)

  # 2. Assess significance
    in_df[,ncol(in_df)+1] = "insign"
    label_sign = paste("gene", sprintf("%03d", which(legend_text < alpha)), sep="")
    in_df[which(in_df[,2]==label_sign),ncol(in_df)] = "signif"
    colnames(in_df)[ncol(in_df)] = "significance"

  # Special to avoid error 'no visible binding for global variable' during compilation
  # For details, see: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    KF_dist = gene_id = significance = NULL

  # 3. Generate plot
    plot_handle = ggplot2::ggplot(data=in_df) +
    geom_density(aes(x=KF_dist, group=gene_id, color=factor(significance), line=2)) +
    xlim(0, xlim_thres) +
    #ylim(0, 25) +
    theme_bw() +
    scale_colour_manual(name="significance", values=rev(color_specs)) +
    ggtitle(paste(title_str, ", alpha=", alpha, "\n", sep="")) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5)))

  # 4. Add annotations
  #    Note: x-position easy to define, because KF distances between 0 and 1
    n_entries = length(legend_text)
    plot_handle = plot_handle + annotate("text", x=annot_x_pos, y=annot_y_pos, label=paste("gene", colnames(legend_text), sep="     "), color=color_specs[2], fontface="bold")
    for (i in 1:n_entries) {
        if (legend_text[i] < alpha) {
            plot_handle = plot_handle + annotate("text", x=annot_x_pos, y=annot_y_pos-(annot_y_pos*i/n_entries), label=paste(rownames(legend_text)[i], sprintf("%03f", legend_text[i]), sep="   "), color=color_specs[1], size=4)
        }
        if (legend_text[i] >= alpha) {
            plot_handle = plot_handle + annotate("text", x=annot_x_pos, y=annot_y_pos-(annot_y_pos*i/n_entries), label=paste(rownames(legend_text)[i], sprintf("%03f", legend_text[i]), sep="   "), color=color_specs[2], size=4)
        }
    }
  # 3. Return plot
    return(plot_handle)
}
