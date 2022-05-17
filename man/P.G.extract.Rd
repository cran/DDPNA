\name{P.G.extract}
\alias{P.G.extract}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Protein Groups information extract.
}
\description{
uniprot ID, ENTRYNAME and status information extract.(only fit for 'MaxQuant' data.)
}
\usage{
P.G.extract(inf, ncol = 4, justID = FALSE,
            status1 = TRUE, ENTRY1 = TRUE, ID1 = TRUE,
            sp1 = TRUE, onlysp = FALSE, verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inf}{
  protein groups IDs information.
}
  \item{ncol}{
  column numbers of output result.
}
  \item{justID}{
  a logical value indicated whether only extract uniprot ID.
}
  \item{status1}{
  a logical value indicated whether extract the first ID status.
}
  \item{ENTRY1}{
  a logical value indicated whether extract the first ID ENTRY NAME.
}
  \item{ID1}{
  a logical value indicated whether extract the first ID UNIPROT ACCESSION.
}
  \item{sp1}{
  a logical value indicated whether extract the first sp ID as the new first ID.
}
  \item{onlysp}{
  a logical value indicated whether only keep sp ID in the second or later ID result.
}
  \item{verbose}{
  integer level of verbosity. Zero means silent, 1 means have Diagnostic Messages.
}
}


\author{
Kefu Liu
}


\examples{
data(ProteomicData)
MaxQdata <- MaxQprotein(ProteomicData$MaxQ)
inf <- P.G.extract(MaxQdata$protein_IDs, justID = TRUE, status = TRUE, ENTRY = TRUE)

}

