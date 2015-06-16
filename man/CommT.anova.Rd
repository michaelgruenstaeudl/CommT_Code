\name{CommT.anova}
\alias{CommT.anova}
\title{Conduct ANOVA analysis on tree distance matrix}

\description{
  Function to conduct an ANOVA analysis on the matrix of tree distances generated via function CommT.kfdist().
}

\usage{
CommT.anova(in_df)
}

\arguments{
  \item{in_df}{a dataframe storing the tree distances inferred via function CommT.kfdist(). This dataframe constitutes a stacked list of distances grouped by genetic locus and is the default output of function CommT.kfdist()}
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
}
