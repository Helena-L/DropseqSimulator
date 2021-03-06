% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_sample_prob_naive.R
\name{calculate_sample_prob_naive}
\alias{calculate_sample_prob_naive}
\title{given a transcript and fragment length, calculate probability of generating reads along the transcript}
\usage{
calculate_sample_prob_naive(transcript_len, fraglen)
}
\arguments{
\item{fraglen}{fragment length, integer.}

\item{transcript}{transcript sequence, DNAString.}
}
\value{
probability of generating reads along the transcript
}
\description{
given a transcript and fragment length, calculate probability of generating reads along the transcript
}
