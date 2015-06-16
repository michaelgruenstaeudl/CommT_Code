\name{CommT-package}
\alias{CommT-package}
\alias{CommT}
\docType{package}

\title{
  Comparative phylogeographic analysis using the community tree framework
}

\description{
  \pkg{CommT} Provides functions to measure the difference between constrained and unconstrained gene tree distributions using various tree distance metrics. Constraints are enforced prior to this analysis via the estimation of a tree under the community tree model.
}

\author{
  Michael Gruenstaeudl

  Maintainer: Michael Gruenstaeudl \email{mi.gruenstaeudl@gmail.com}
}

\references{
  Gruenstaeudl, M., Reid, N.M., Wheeler, G.R. and Carstens, B.C., 2015. Posterior Predictive Checks of Coalescent Models: P2C2M, an R package. Molecular Ecology Resources, in press.
}

\examples{
  # Load libraries
  #  require(ape)

  # Load input data
    data(post_gt_distrs_BEAST)
    data(post_gt_distrs_starBEAST)

  # Calculate KF distances
    in_data = CommT.kfdist(post_gt_distrs_BEAST, post_gt_distrs_starBEAST)

  # Generate ANOVA legend
    legend_text = CommT.anova(in_data)

  # Generate coordinates for plot
    legend_pos = CommT.legendpos(in_data)

  # Visualize KF distances
    my_plot = CommT.viz(in_data, "my_project", alpha=0.05, legend_text, legend_pos)

  # Display plot
    my_plot
}
