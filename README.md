***Note:* For questions and discussion about [msisensor](https://github.com/ding-lab/msisensor), please visit the repository of [msisensor-pro](https://github.com/xjtu-omics/msisensor-pro) at https://github.com/xjtu-omics/msisensor-pro.**



## MSIsensor



**MSIsensor** is a C++ program to detect replication slippage variants at microsatellite regions, and differentiate them as somatic or germline. Given paired tumor and normal sequence data, it builds a distribution for expected (normal) and observed (tumor) lengths of repeated sequence per microsatellite, and compares them using Pearson's Chi-Squared Test. Comprehensive testing indicates MSIsensor is an efficient and effective tool for deriving microsatellite instability (MSI) status from standard tumor-normal paired sequence data. MSIsensor is publiched in [*Bioinformatics*](https://www.ncbi.nlm.nih.gov/pubmed/24371154). Please click [here](https://github.com/PengJia6/msisensor/blob/master/README_MSIsensor.md) to see more details about MSIsensor. If you have any questions about MSIsensor, please contact one or more of the following folks:
Beifang Niu (<bniu@sccas.cn>), Kai Ye (<kaiye@xjtu.edu.cn>) or Li Ding (<lding@wustl.edu>).

If you used this tool for your work, please cite [PMID 24371154](https://www.ncbi.nlm.nih.gov/pubmed/24371154)


**[MSIsensor-pro](https://github.com/xjtu-omics/msisensor-pro)** is a new MSI detection method developed by [Kai Ye](http://gr.xjtu.edu.cn/web/kaiye/english)  et al. [MSIsensor-pro](https://github.com/xjtu-omics/msisensor-pro) is a fast, accurate, and matched-normal-sample-free MSI detection method. It accepts the whole genome sequencing, whole exome sequencing and target region (panel) sequencing data as input. MSIsensor-pro introduces a multinomial distribution model to quantify polymerase slippages for each tumor sample and a discriminative sites selection method to enable MSI detection without matched normal samples. **[MSIsensor-pro](https://github.com/xjtu-omics/msisensor-pro)** is now published in [*Genomics Proteomics & Bioinformatics*](https://www.sciencedirect.com/science/article/pii/S1672022920300218). If you have any question about MSIsensor-pro, please [open a issue](https://github.com/xjtu-omics/msisensor-pro/issues/new) on [MSIsensor-pro's homepage](https://github.com/xjtu-omics/msisensor-pro) or contact with Kai Ye (kaiye@xjtu.edu.cn) directly. 


**MSIsensor2** is also a MSI detecteion method specially designed for tumor only sequencing data. MSIsensor2 was developed by Beifang Niu's lab (bniu@sccas.cn) independently. Please try the MSIsensor2 here: https://github.com/niu-lab/msisensor2 or require any further details here: http://niulab.scgrid.cn/msisensor2/index.html.

