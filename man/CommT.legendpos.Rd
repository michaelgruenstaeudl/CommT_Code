\name{CommT.legendpos}
\alias{CommT.legendpos}
\title{Infer coordinates for plot legend}

\description{
  Function to infer the optimal position of the ANOVA legend given the calculated tree distance distributions.
}

\usage{
CommT.legendpos(in_data)
}

\arguments{
  \item{in_data}{a dataframe storing the tree distances inferred via function CommT.kfdist(). This dataframe constitutes a stacked list of distances grouped by genetic locus and is the default output of function CommT.kfdist()}
}

\author{
  Michael Gruenstaeudl

  Maintainer: Michael Gruenstaeudl \email{mi.gruenstaeudl@gmail.com}
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

}
