\name{CommT.kfdist}
\alias{CommT.kfdist}
\title{calculate the tree distances}

\description{
  This function calculates the tree distances.
}

\usage{
CommT.kfdist(post_gt_distrs_BEAST, post_gt_distrs_starBEAST,
             outlier_num = 1, treedist_select = 2)
}

\arguments{
  \item{post_gt_distrs_BEAST}{Bla bla bla.}

  \item{post_gt_distrs_starBEAST}{Bla bla bla.}

  \item{outlier_num}{an integer. The default is \code{1}.}

  \item{treedist_select}{an integer. The default is \code{2}.}
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
}
