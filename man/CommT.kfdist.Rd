\name{CommT.kfdist}
\alias{CommT.kfdist}
\title{Calculate the tree distance metric}

\description{
  Function to calculate a selected tree distance metric between posterior gene tree distributions.
}

\usage{
CommT.kfdist(post_gt_distrs_BEAST, post_gt_distrs_starBEAST,
             outlier_num = 1, treedist_select = 2)
}

\arguments{
  \item{post_gt_distrs_BEAST}{a multiPhylo object representing the posterior gene tree distribution estimated in the \bold{absence} of a multispecies coalescent model. This object is split into individual lists of trees, whereby each list represents a posterior gene tree distribution inferred under a specific genetic locus. Each list should be named by the respective locus.}

  \item{post_gt_distrs_starBEAST}{a multiPhylo object representing the posterior gene tree distribution estimated in the \bold{presence} of a multispecies coalescent model. This object is split into individual lists of trees, whereby each list represents a posterior gene tree distribution inferred under a specific genetic locus. Each list should be named by the respective locus.}

  \item{outlier_num}{an integer specifying the number of the genetic locus to be used as positive control in the ANOVA calculations employing mixed effects models. The default is \code{1}.}

  \item{treedist_select}{an integer specifying the type of tree distance metric. The selection hereby follows the output order employed in function \code{treedist} of the R package \code{phangorn}. The default is \code{2}, which represents the Kuhner-Felsenstein distance.}
}

\author{
  Michael Gruenstaeudl

  Maintainer: Michael Gruenstaeudl \email{mi.gruenstaeudl@gmail.com}
}

\references{
  Kuhner, M.K. and Felsenstein. J., 1994. Simulation comparison of phylogeny algorithms under equal and unequal evolutionary rates. Molecular Biology and Evolution 11: 459-468.
}

\examples{
  # Load libraries
  #  require(ape)

  # Load input data
    data(post_gt_distrs_BEAST)
    data(post_gt_distrs_starBEAST)
    
  # Calculate KF distances
    in_data = CommT.kfdist(post_gt_distrs_BEAST, post_gt_distrs_starBEAST)
}
