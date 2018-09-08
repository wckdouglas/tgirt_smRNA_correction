# Model set up #

We hypothesized the first and last three bases of the RNA fragments are the sources of bias in TGIRT-seq. We assume that each positional nucleotide has a independent weight ($p_{i,b}$), where $i\in(1,2,3)$ and $b\in{A,C,T,G}$. We can model the observed count-pre-million ($CPM_{obs}$) values with true CPM ($CPM_{true}$) values by:

$$ CPM_{obs} = CPM_{true} \prod_{i=1,2,3}\prod_{b\in(A,C,T,G)}p_{i,b} $$

$$ \frac{CPM_{obs}}{CPM_{true}} = \prod_{i}\prod_{b}p_{i,b} $$ 

$$ log(\frac{CPM_{obs}}{CPM_{true}}) = \sum_{i}\sum_{b}log(p_{i,b}) $$ 

$$ log(CPM_{obs}) - log(CPM_{true}) = \sum_{i}\sum_{b}log(p_{i,b})  $$

The weights can be derived by Ridge regression using one-hot-encoded nucleotide sequences in the form of:

$$ \Delta log (CPM) = \sum_{i}\sum_{b}log(p_{i,b}) $$

Suppose we have a miRNA starting with 5'-AGG,

$$ \Delta log (CPM) = (1)log(p_{1,A}) + (0)log(p_{1,G}) + (0)log(p_{1,T}) +  ... $$



# Correction #

We define read count (RC) values as :

$$ CPM = \frac{10^6 RC}{\sum RC} $$

Thus, 
$$ RC_{true} = 10^6  CPM_{true} \sum RC_{true}$$
$$ RC_{obs} = 10^6  CPM_{obs} \sum RC_{obs}$$

To get the true read count ($RC_{true}$)

$$ log(CPM_{obs}) - log(CPM_{true}) = \sum_{i}\sum_{b}log(p_{i,b})  $$
$$ log(CPM_{true}) = log(CPM_{obs}) - \sum_{i}\sum_{b}log(p_{i,b}) $$
$$ log(\frac{10^6 RC_{true}}{\sum RC_{true}}) = log(\frac{10^6 RC_{obs}}{\sum RC_{obs}}) - \sum_{i}\sum_{b}log(p_{i,b}) $$
$$ log(\frac{RC_{true}}{\sum RC_{true}}) + log(10^6) = log(\frac{RC_{obs}}{\sum RC_{obs}}) + log({10^6}) - \sum_{i}\sum_{b}log(p_{i,b}) $$
$$ log(\frac{RC_{true}}{\sum RC_{true}}) = log(\frac{RC_{obs}}{\sum RC_{obs}}) - \sum_{i}\sum_{b}log(p_{i,b}) $$
$$ \frac{RC_{true}}{\sum RC_{true}} = exp^{(log(\frac{RC_{obs}}{\sum RC_{obs}}) - \sum_{i}\sum_{b}log(p_{i,b}))} $$

