\name{exact_zeroindex}
\alias{exact_zeroindex}

\title{
Exact indices of zero elements
}
\description{
This function exact indices of entries of support matrix to be constrained to be zero.
}
\usage{
exact_zeroindex(support = NULL)
}
\arguments{
  \item{support}{
  A p*p supprt matrix. The zero elements within the support matrix represent the disconnected edges in GRN.
}
}

\value{
  A two-column matrix, and each row indicating the indices of elements in GRN to be constrained to be zero. Furthermore, if the input "support" is null, the output is null.
}

\examples{
##1. Package loading
##-------------------------------------------------------------
library(scGeneNet)
##-------------------------------------------------------------

##2. Generate the used support
##-------------------------------------------------------------
support_use<-diag(100)
support_use[1:10,]<-1
support_use[,1:10]<-1
colnames(support_use)<-paste("Gene",1:100,sep = "")
rownames(support_use)<-paste("Gene",1:100,sep = "")
##-------------------------------------------------------------

##3. Induce the index of zero element
##-------------------------------------------------------------
## (a) give used support
zero_GRN<-exact_zeroindex(support = support_use)

## (b) not give used support
zero_GRN_null<-exact_zeroindex(support = NULL)
##-------------------------------------------------------------
}