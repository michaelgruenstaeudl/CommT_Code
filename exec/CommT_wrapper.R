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

      # Keywords: noswap, g1swap, g1a2swap, g1t4swap, gAllShuf

      # Specify input file names (User input necessary)
        BEAST_prefix = paste("/home/michael_science/Desktop/CommT_Project/Revision_SupplFigure1/05_Project_CommT_Jan2015/03_KFdist_distributions/01_input/01_starBEAST_geneTrees/CommT.Jan2015.g1t4swap.sim.", s, ".gene", sep="")
        starBEAST_prefix = paste("/home/michael_science/Desktop/CommT_Project/Revision_SupplFigure1/05_Project_CommT_Jan2015/03_KFdist_distributions/01_input/02_BEAST_GeneTrees/CommT.Jan2015.g1t4swap.sim.", s, ".gene", sep="")

      # Load data
        post_gt_distrs_BEAST = list()
        post_gt_distrs_starBEAST = list()
        for (g in 1:n_genes) {
            g_lz = sprintf("%03d", g)
            post_gt_distrs_BEAST[[g]] = read.nexus(paste(BEAST_prefix, g_lz, ".trees", sep=""))
            post_gt_distrs_starBEAST[[g]] = read.nexus(paste(starBEAST_prefix, g_lz, ".trees", sep=""))
        }

      # Calculate KF distances
        in_data = CommT.kfdist(post_gt_distrs_BEAST, post_gt_distrs_starBEAST, outlier_num=c())

      # Generate ANOVA legend
        legend_text = CommT.anova(in_data)

      # Generate coordinates for plot
        legend_pos = CommT.legendpos(in_data)

      # Visualize KF distances
        #my_plot = CommT.viz(in_data, "a_project_name_here", alpha=0.05, legend_text, legend_pos)
        my_plot = CommT.viz(in_data, paste("sim.", s, sep=""), alpha=0.05, legend_text, legend_pos)

      # Save plot
        assign(paste("sim.", s, sep=""), my_plot)

    }

    # Specify output file names (User input necessary)
    svg("/home/michael_science/Desktop/CommT_Project/Revision_SupplFigure1/05_Project_CommT_Jan2015/03_KFdist_distributions/03_original_visualizations_with_Rpackage/05_BEAST.to.starBEAST_KFdist_g1t4swap.svg", width=30, height=26)

    grid.arrange(sim.01, sim.02, sim.03, sim.04, sim.05,
                 sim.06, sim.07, sim.08, sim.09, sim.10,
                 sim.11, sim.12, sim.13, sim.14, sim.15,
                 sim.16, sim.17, sim.18, sim.19, sim.20, ncol=5)

    dev.off()
