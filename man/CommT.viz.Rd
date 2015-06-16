\name{CommT.viz}
\alias{CommT.viz}
\title{Visualize tree distances and ANOVA results}

\description{
  Function to generate a visualization of both the tree distances and ANOVA results.
}

\usage{
CommT.viz(in_df, title_str="a_project_name_here", alpha=0.05, legend_text, legend_pos)
}

\arguments{
  \item{in_df}{a dataframe storing the tree distances inferred via function CommT.kfdist(). This dataframe constitutes a stacked list of distances grouped by genetic locus and is the default output of function CommT.kfdist()}

  \item{title_str}{a string specifying a title for the resulting plot. The default is \code{a_project_name_here}.}

  \item{alpha}{a float specifying the alpha-level employed in the ANOVA. The default is \code{0.05}.}

  \item{legend_text}{a list specifying the ANOVA results to plotted. It constitutes the default output of function CommT.anova()}
  
  \item{legend_pos}{a list specifying the plot position. It constitutes the default output of function CommT.legendpos()}

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

  # Visualize KF distances
    my_plot = CommT.viz(in_data, "my_project", alpha=0.05, legend_text, legend_pos)
}
