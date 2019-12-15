# Repository for manuscript "Autosomal genetic variation is associated with DNA methylation in regions variably escaping X-chromosome inactivation" by Luijk _et al_.

#### Calculate P-values for all SNP-CpG pairs

_MatrixEQTL.R_ loads the data and fits a linear regression model for every SNP-CpG pair using MatrixEQTL. The variable _gender_ has to be set to either 'females' or 'males'. The output will be written per gender.

#### Calculate overall Simes P-value per SNP

_Simes.R_ read the output from MatrixEQTL and uses the calculated P-values per SNP-CpG pair to calculate one overall P-value per SNP. Again, the variable _gender_ has to be set.
