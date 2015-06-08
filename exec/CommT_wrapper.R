    library(ape)
    library(CommT)
    library(ggplot2)
    library(gridExtra)
    library(phangorn)
    library(reshape)

    n_genes = 10
    n_sims = sprintf("%02d", 1:20)

    for (s in n_sims) {
        cat(paste("\n", "Analyzing simulation", s, "\n"))

      # Set prefices
      
        #BEAST_prefix = paste("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/01_input/02_BEAST_GeneTrees/CommT.Jan2015.noswap.sim.", s, ".gene", sep="")
        #starBEAST_prefix = paste("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/01_input/01_starBEAST_geneTrees/CommT.Jan2015.noswap.sim.", s, ".gene", sep="")
        
        #BEAST_prefix = paste("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/01_input/02_BEAST_GeneTrees/CommT.Jan2015.g1swap.sim.", s, ".gene", sep="")
        #starBEAST_prefix = paste("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/01_input/01_starBEAST_geneTrees/CommT.Jan2015.g1swap.sim.", s, ".gene", sep="")
        
        #BEAST_prefix = paste("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/01_input/02_BEAST_GeneTrees/CommT.Jan2015.g1a2swap.sim.", s, ".gene", sep="")
        #starBEAST_prefix = paste("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/01_input/01_starBEAST_geneTrees/CommT.Jan2015.g1a2swap.sim.", s, ".gene", sep="")

        #BEAST_prefix = paste("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/01_input/02_BEAST_GeneTrees/CommT.Jan2015.g1t4swap.sim.", s, ".gene", sep="")
        #starBEAST_prefix = paste("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/01_input/01_starBEAST_geneTrees/CommT.Jan2015.g1t4swap.sim.", s, ".gene", sep="")

        BEAST_prefix = paste("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/01_input/05_gAllShuf/BEAST/CommT.Jan2015.gAllShuf.sim.", s, ".gene", sep="")
        starBEAST_prefix = paste("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/01_input/05_gAllShuf/starBEAST/CommT.Jan2015.gAllShuf.sim.", s, ".gene", sep="")

      # Load data
        post_gt_distrs_BEAST = list()
        post_gt_distrs_starBEAST = list()
        for (g in 1:n_genes) {
            g_lz = sprintf("%03d", g)
            post_gt_distrs_BEAST[[g]] = read.nexus(paste(BEAST_prefix, g_lz, ".trees", sep=""))
            post_gt_distrs_starBEAST[[g]] = read.nexus(paste(starBEAST_prefix, g_lz, ".trees", sep=""))
        }

      # Calculate KF distances
        #in_data = CommT.kfdist(post_gt_distrs_BEAST, post_gt_distrs_starBEAST, outlier_num=c())
        #in_data = CommT.kfdist(post_gt_distrs_BEAST, post_gt_distrs_starBEAST, outlier_num=1)
        #in_data = CommT.kfdist(post_gt_distrs_BEAST, post_gt_distrs_starBEAST, outlier_num=c(1,2))
        #in_data = CommT.kfdist(post_gt_distrs_BEAST, post_gt_distrs_starBEAST, outlier_num=c(1,2,3,4))
        in_data = CommT.kfdist(post_gt_distrs_BEAST, post_gt_distrs_starBEAST, outlier_num=c())
        
      # Generate ANOVA legend
        in_legend = CommT.anova(in_data)
    
      # Correct legend
        rownames(in_legend) = post_fxs
        
      # Generate coordinates for plot
        annot_x_list = annot_y_list = xlim_thres_list = c()
        for (i in split(in_data, in_data$gene_id)) {
            xlim_thres_list = c(annot_x_list, quantile(i[,"KF_dist"], probs = c(0.9999)))
            annot_y_list = c(annot_y_list, max(table(cut(i[,"KF_dist"], breaks=100))))
        }
        xlim_thres_pos = quantile(xlim_thres_list, probs = c(0.95))
        annot_x_pos = quantile(xlim_thres_list, probs = c(0.80))[[1]]
        annot_y_pos = max(annot_y_list)

      # Visualize KF distances
        my_plot = CommT.viz(in_data, in_legend, title_str=paste("sim.", s, sep=""),
                            alpha=0.05, annot_x_pos, annot_y_pos, xlim_thres_pos)

      # Save plot
        assign(paste("sim.", s, sep=""), my_plot)

    }

    #svg("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/02_visualizations_with_Rpackage/01_BEAST.to.starBEAST_KFdist_noswap.svg", width=30, height=20)
    #svg("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/02_visualizations_with_Rpackage/02_BEAST.to.starBEAST_KFdist_g1swap.svg", width=30, height=20)
    #svg("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/02_visualizations_with_Rpackage/03_BEAST.to.starBEAST_KFdist_g1a2swap.svg", width=30, height=20)
    #svg("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/02_visualizations_with_Rpackage/04_BEAST.to.starBEAST_KFdist_g1t4swap.svg", width=30, height=20)
    svg("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/02_visualizations_with_Rpackage/05_BEAST.to.starBEAST_KFdist_gAllShuf.svg", width=30, height=20)

    grid.arrange(sim.01, sim.02, sim.03, sim.04, sim.05,
                 sim.06, sim.07, sim.08, sim.09, sim.10,
                 sim.11, sim.12, sim.13, sim.14, sim.15,
                 sim.16, sim.17, sim.18, sim.19, sim.20, ncol=5)
    dev.off()
