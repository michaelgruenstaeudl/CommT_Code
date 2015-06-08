\name{CommT-package}
\alias{CommT-package}
\alias{CommT}
\docType{package}

\title{
  Evaluations of the community tree
}

\description{
  \pkg{CommT} Placeholder Placeholder Placeholder.
}

\note{

  \bold{Placeholder}

  To use \pkg{CommT}, Placeholder Placeholder Placeholder:

  \code{echo 'alias python=python2.7' >> ~/.bashrc}
}

\author{
  Michael Gruenstaeudl

  Maintainer: Michael Gruenstaeudl \email{gruenstaeudl.1@osu.edu}
}

\references{
  Gruenstaeudl, M., Reid, N.M., Wheeler, G.R. and Carstens, B.C., submitted. Posterior Predictive Checks of Coalescent Models: P2C2M, an R package.
}

\examples{
  # Load input data
    data(post_gt_distrs_BEAST)
    data(post_gt_distrs_starBEAST)

  # Calculate KF distances
    in_data = CommT.kfdist(post_gt_distrs_BEAST, post_gt_distrs_starBEAST, outlier_num=1)

  # Generate ANOVA legend
    in_legend = CommT.anova(in_data)

  # Visualize KF distances
    my_plot = CommT.viz(in_data, in_legend, "test simulation")
}
