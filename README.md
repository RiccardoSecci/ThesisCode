# ThesisCode
In this repository, I put all the basis code used in order to obtain the results of my Masters Thesis.
It also contains the Results of the Code. From the Final Heatmap that was used in order to identify Drugs that 
promote Hematopoietic Differentiation, to all the side results used to further understand and interpret the Results.

Basically, my Project start from the Connectivity Map Database, and uses a Transcription Factor tool to Find the TF that are 
enriched by each Compound. Then, the user can select a TF Subset, in my case I chose to focus on TFs related to Hematopoietic Differentiation at first, to identify pro-Hematopoietic drugs for the treatment of Acute Myeloid Leukemia. 
I then did a second Analysis using TFs related to Apoptosis, in order to select out of the top scoring drugs calculated in the previous round the one that also caused Apoptosis.
That's because I (and my team) was looking for Drugs that could be suitable for a Pro-Differentiative therapy, not a Cytotoxic Chemotherapy.
I also did a lot of sideline Analyses, which I explained thoroughly in my thesis, from dimensionality reduction methods to understand the gene expression data to Enrichment for Kegg, Go and Wikipathways terms, in order to have other instruments for the Analysis of my Results.

In the end, the top scoring drugs are being analysing by a PhD student from the lab I worked in, and the early results look very promising.
Known drugs that have a pro-differentiative effect appear in the top scoring ones, and the new identified pro-differentiative drugs seem to have this property, at leat in vitro.

I coded both in Python and in R, as some tasks were easier with the former, some with the latter.

