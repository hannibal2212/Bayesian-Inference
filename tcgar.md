``` r
library(readr)
```

    ## Warning: package 'readr' was built under R version 4.2.2

``` r
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 4.2.3

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(TCGAbiolinks)
library(stringr)
```

    ## Warning: package 'stringr' was built under R version 4.2.2

``` r
library(magrittr)
```

    ## Warning: package 'magrittr' was built under R version 4.2.1

## primary query for looking at the number of samples in the TCGA database

``` r
primary_query<-GDCquery(project = "TCGA-BRCA",data.category = "Transcriptome Profiling",
         data.type = "Gene Expression Quantification",
         workflow.type = "STAR - Counts",
         sample.type = c("Primary Tumor", "Solid Tissue Normal")
         )
```

    ## --------------------------------------

    ## o GDCquery: Searching in GDC database

    ## --------------------------------------

    ## Genome of reference: hg38

    ## --------------------------------------------

    ## oo Accessing GDC. This might take a while...

    ## --------------------------------------------

    ## ooo Project: TCGA-BRCA

    ## --------------------

    ## oo Filtering results

    ## --------------------

    ## ooo By data.type

    ## ooo By workflow.type

    ## ooo By sample.type

    ## ----------------

    ## oo Checking data

    ## ----------------

    ## ooo Checking if there are duplicated cases

    ## ooo Checking if there are results for the query

    ## -------------------

    ## o Preparing output

    ## -------------------

``` r
# number of cases and controls totally availabe in the TCGA
primary_query[[1]][[1]]%>%group_by(sample_type)%>%summarise(n = n())
```

    ## # A tibble: 2 × 2
    ##   sample_type             n
    ##   <chr>               <int>
    ## 1 Primary Tumor        1111
    ## 2 Solid Tissue Normal   113

## reading the sample data file

``` r
sample_data<-read.csv("sampledata.csv", stringsAsFactors = F)
```

``` r
head(sample_data)
```

    ##   X     TCGA_SAMPLE      BARCODE PAM50 PAM50lite TNBCtype TNBCtype_4 gender
    ## 1 1 TCGA-A1-A0SP-01 TCGA-A1-A0SP Basal     Basal      UNC        BL1 female
    ## 2 2 TCGA-A2-A0YM-01 TCGA-A2-A0YM Basal     Basal      BL1        BL1 female
    ## 3 3 TCGA-A2-A3XS-01 TCGA-A2-A3XS Basal     Basal      BL1        BL1 female
    ## 4 4 TCGA-A2-A3XT-01 TCGA-A2-A3XT Basal     Basal      BL1        BL1 female
    ## 5 5 TCGA-A2-A3XX-01 TCGA-A2-A3XX Basal     Basal      BL1        BL1 female
    ## 6 6 TCGA-AO-A124-01 TCGA-AO-A124 Basal     Basal      BL1        BL1 female
    ##   primarysiteofdisease neoplasm.diseasestage pathology.T.stage
    ## 1               breast             stage iia                t2
    ## 2               breast             stage iia                t2
    ## 3               breast            stage iiia                t1
    ## 4               breast             stage iib                t2
    ## 5               breast             stage iia                t2
    ## 6               breast             stage iia                t2
    ##   pathology.N.stage pathology.M.stage number.of.lymph.nodes vitalstatus_TNBC
    ## 1           n0 (i-)                m0                     0                0
    ## 2           n0 (i-)                m0                     0                0
    ## 3               n2a                m0                     9                1
    ## 4               n1a                m0                     1                0
    ## 5                n0                m0                     0                1
    ## 6           n0 (i-)                m0                     0                0
    ##   yearstodeath_TNBC              histologicaltype ESTIMATE ABSOLUTE   LUMP
    ## 1          1.600000 infiltrating ductal carcinoma   0.8060     0.46 0.6603
    ## 2          2.643836 infiltrating ductal carcinoma   0.7678     0.58 0.7495
    ## 3          2.827397 infiltrating ductal carcinoma   0.8091       NA 0.7879
    ## 4          6.917808 infiltrating ductal carcinoma   0.6980       NA 0.6838
    ## 5          3.942466 infiltrating ductal carcinoma   0.8138       NA 0.7889
    ## 6          8.547945                other  specify   0.9082     0.81 0.6698
    ##     IHC    CPE CHAT.purity NMF..Gene. NMF..miRNA. NMF..CNV. SNF..gene...miRNA.
    ## 1 0.725 0.6913        0.58          1           1         1                  1
    ## 2 0.910 0.7266        0.76          3           2         2                  1
    ## 3 0.680 0.7855        0.69          2           4         3                  1
    ## 4 0.850 0.6632        0.57          3           1         2                  1
    ## 5 0.880 0.7868          NA          3           2         3                  1
    ## 6 0.725 0.8859        0.87          2           3         3                  1
    ##   SNF..gene...miRNA...CNV. SNF..gene...miRNA...CNV...meth.
    ## 1                        3                               1
    ## 2                        1                               1
    ## 3                        1                               1
    ## 4                        1                               1
    ## 5                        1                               1
    ## 6                        1                               1
    ##                                     id data_format                        cases
    ## 1 222b562a-3ef7-4e69-a60d-54aac5104102         TSV TCGA-A1-A0SP-01A-11R-A084-07
    ## 2 fe7dfda5-6846-4238-9ccc-472978eb78a1         TSV TCGA-A2-A0YM-01A-11R-A109-07
    ## 3 14979d31-e8b1-4948-beea-f913dfefa200         TSV TCGA-A2-A3XS-01A-11R-A22U-07
    ## 4 0c068d5a-ccad-47a5-bfe0-0dc5f136d6aa         TSV TCGA-A2-A3XT-01A-11R-A22U-07
    ## 5 c0f6d246-c872-41a5-8135-2bdbdb1f924f         TSV TCGA-A2-A3XX-01A-21R-A239-07
    ## 6 a63cca36-139d-4054-9d99-9b04ddae7ec2         TSV TCGA-AO-A124-01A-11R-A10J-07
    ##   access
    ## 1   open
    ## 2   open
    ## 3   open
    ## 4   open
    ## 5   open
    ## 6   open
    ##                                                                     file_name
    ## 1 e0527b4b-de23-4b9f-8871-13f64f6976cd.rna_seq.augmented_star_gene_counts.tsv
    ## 2 8a27980e-a506-4cb3-91f2-3f0e5a19acfe.rna_seq.augmented_star_gene_counts.tsv
    ## 3 deb14c1d-474c-4d6f-bb8a-2b6df5dce5c9.rna_seq.augmented_star_gene_counts.tsv
    ## 4 bdba14fb-88aa-43b0-8c25-53948b9935b7.rna_seq.augmented_star_gene_counts.tsv
    ## 5 f9506da9-20a3-4124-a614-eb525787ea45.rna_seq.augmented_star_gene_counts.tsv
    ## 6 3d05d072-d042-4c02-aaf6-ee2e589003d6.rna_seq.augmented_star_gene_counts.tsv
    ##                           submitter_id           data_category            type
    ## 1 01c7a57e-080f-4e18-ba49-c80242b09a33 Transcriptome Profiling gene_expression
    ## 2 d9f79023-a752-4d31-abc7-694221916611 Transcriptome Profiling gene_expression
    ## 3 2596a9d3-894b-4466-80ed-3181500b9a5c Transcriptome Profiling gene_expression
    ## 4 6531f2c8-bdb7-423e-bb35-7f1b64d4717e Transcriptome Profiling gene_expression
    ## 5 ead405b4-03cb-43c8-bb1b-aada1eb5cebc Transcriptome Profiling gene_expression
    ## 6 2db4cd6c-30c8-4c91-afe6-9b3ad45710b0 Transcriptome Profiling gene_expression
    ##   file_size                 created_datetime                           md5sum
    ## 1   4244140 2021-12-13T22:43:37.520668-06:00 d4c737d27e18342109001a0e7981dcfe
    ## 2   4249051 2021-12-13T22:06:39.713399-06:00 5ab2e6521e4cfa4a63dec5b21bc66ff4
    ## 3   4252274 2021-12-13T22:08:40.155876-06:00 36946ef0c730515335ad58f52942f433
    ## 4   4249842 2021-12-13T21:46:40.906869-06:00 518eace9d4d12dc66d378770abbb0f00
    ## 5   4248329 2021-12-13T22:20:08.151779-06:00 03161539d6622a38903b6e76b3ce92c0
    ## 6   4249197 2021-12-13T21:55:15.782366-06:00 f2f1a4fff2ed21eac2b8fb473d13a10b
    ##                   updated_datetime                              file_id
    ## 1 2022-01-19T12:00:30.311996-06:00 222b562a-3ef7-4e69-a60d-54aac5104102
    ## 2 2022-01-19T12:00:30.311996-06:00 fe7dfda5-6846-4238-9ccc-472978eb78a1
    ## 3 2022-01-19T12:00:30.311996-06:00 14979d31-e8b1-4948-beea-f913dfefa200
    ## 4 2022-01-19T12:00:30.311996-06:00 0c068d5a-ccad-47a5-bfe0-0dc5f136d6aa
    ## 5 2022-01-19T12:00:30.311996-06:00 c0f6d246-c872-41a5-8135-2bdbdb1f924f
    ## 6 2022-01-19T12:00:30.311996-06:00 a63cca36-139d-4054-9d99-9b04ddae7ec2
    ##                        data_type    state experimental_strategy version
    ## 1 Gene Expression Quantification released               RNA-Seq       1
    ## 2 Gene Expression Quantification released               RNA-Seq       1
    ## 3 Gene Expression Quantification released               RNA-Seq       1
    ## 4 Gene Expression Quantification released               RNA-Seq       1
    ## 5 Gene Expression Quantification released               RNA-Seq       1
    ## 6 Gene Expression Quantification released               RNA-Seq       1
    ##   data_release   project                          analysis_id analysis_state
    ## 1  32.0 - 38.0 TCGA-BRCA c6b1e0f7-2ca0-427b-a628-7679365e5a76       released
    ## 2  32.0 - 38.0 TCGA-BRCA 8e5789b2-f949-4d46-bee6-34d1907c37f2       released
    ## 3  32.0 - 38.0 TCGA-BRCA 69a76a2b-5923-4955-b228-e794675a4722       released
    ## 4  32.0 - 38.0 TCGA-BRCA 0db3347c-d411-4fe6-99cb-61b51427a821       released
    ## 5  32.0 - 38.0 TCGA-BRCA 0130ef88-c00c-4e68-9337-84c22546600f       released
    ## 6  32.0 - 38.0 TCGA-BRCA f205b26c-ce20-4456-a8b2-b8ba9ca03f47       released
    ##                               analysis_submitter_id
    ## 1 e0527b4b-de23-4b9f-8871-13f64f6976cd_star__counts
    ## 2 8a27980e-a506-4cb3-91f2-3f0e5a19acfe_star__counts
    ## 3 deb14c1d-474c-4d6f-bb8a-2b6df5dce5c9_star__counts
    ## 4 bdba14fb-88aa-43b0-8c25-53948b9935b7_star__counts
    ## 5 f9506da9-20a3-4124-a614-eb525787ea45_star__counts
    ## 6 3d05d072-d042-4c02-aaf6-ee2e589003d6_star__counts
    ##                                                                                                                         analysis_workflow_link
    ## 1 https://github.com/NCI-GDC/gdc-rnaseq-cwl/blob/5d8c131bbff59fb0c969217fc1d44e6d1503cd1f/rnaseq-star-align/star2pass.rnaseq_harmonization.cwl
    ## 2 https://github.com/NCI-GDC/gdc-rnaseq-cwl/blob/5d8c131bbff59fb0c969217fc1d44e6d1503cd1f/rnaseq-star-align/star2pass.rnaseq_harmonization.cwl
    ## 3 https://github.com/NCI-GDC/gdc-rnaseq-cwl/blob/5d8c131bbff59fb0c969217fc1d44e6d1503cd1f/rnaseq-star-align/star2pass.rnaseq_harmonization.cwl
    ## 4 https://github.com/NCI-GDC/gdc-rnaseq-cwl/blob/5d8c131bbff59fb0c969217fc1d44e6d1503cd1f/rnaseq-star-align/star2pass.rnaseq_harmonization.cwl
    ## 5 https://github.com/NCI-GDC/gdc-rnaseq-cwl/blob/5d8c131bbff59fb0c969217fc1d44e6d1503cd1f/rnaseq-star-align/star2pass.rnaseq_harmonization.cwl
    ## 6 https://github.com/NCI-GDC/gdc-rnaseq-cwl/blob/5d8c131bbff59fb0c969217fc1d44e6d1503cd1f/rnaseq-star-align/star2pass.rnaseq_harmonization.cwl
    ##   analysis_workflow_type                analysis_workflow_version   sample_type
    ## 1          STAR - Counts 5d8c131bbff59fb0c969217fc1d44e6d1503cd1f Primary Tumor
    ## 2          STAR - Counts 5d8c131bbff59fb0c969217fc1d44e6d1503cd1f Primary Tumor
    ## 3          STAR - Counts 5d8c131bbff59fb0c969217fc1d44e6d1503cd1f Primary Tumor
    ## 4          STAR - Counts 5d8c131bbff59fb0c969217fc1d44e6d1503cd1f Primary Tumor
    ## 5          STAR - Counts 5d8c131bbff59fb0c969217fc1d44e6d1503cd1f Primary Tumor
    ## 6          STAR - Counts 5d8c131bbff59fb0c969217fc1d44e6d1503cd1f Primary Tumor
    ##   is_ffpe sample.submitter_id
    ## 1   FALSE    TCGA-A1-A0SP-01A
    ## 2   FALSE    TCGA-A2-A0YM-01A
    ## 3   FALSE    TCGA-A2-A3XS-01A
    ## 4   FALSE    TCGA-A2-A3XT-01A
    ## 5   FALSE    TCGA-A2-A3XX-01A
    ## 6   FALSE    TCGA-AO-A124-01A

``` r
head(primary_query[[1]][[1]])
```

    ##                                     id data_format                        cases
    ## 1 a5c799bc-1332-4127-8fda-84da04084d25         TSV TCGA-A7-A26E-01B-06R-A277-07
    ## 2 1dd898ab-3b85-4932-b889-0eb0a494f2b6         TSV TCGA-A2-A0CU-01A-12R-A034-07
    ## 3 7286a685-17ef-4ef5-9b7b-2696c6e7e1a6         TSV TCGA-PL-A8LV-01A-21R-A41B-07
    ## 4 de6f1503-33d7-4f86-b835-bdffba7ea4e3         TSV TCGA-BH-A0BC-01A-22R-A084-07
    ## 5 4c46a463-f0e3-4e19-ad95-5dab0231d1c6         TSV TCGA-AR-A1AX-01A-11R-A12P-07
    ## 6 586d9c21-e58a-4171-9755-30c0e8b16293         TSV TCGA-AC-A2FO-01A-11R-A180-07
    ##   access
    ## 1   open
    ## 2   open
    ## 3   open
    ## 4   open
    ## 5   open
    ## 6   open
    ##                                                                     file_name
    ## 1 de01516e-43f0-4f96-8ac6-ab543a314829.rna_seq.augmented_star_gene_counts.tsv
    ## 2 e0e050cc-5e2a-41fc-8394-b4816a9dc8a2.rna_seq.augmented_star_gene_counts.tsv
    ## 3 3572211c-0611-4d67-9e71-c5a8f704bd1f.rna_seq.augmented_star_gene_counts.tsv
    ## 4 195d8bd4-352c-4407-9319-797d915c72b8.rna_seq.augmented_star_gene_counts.tsv
    ## 5 f2ab4258-222d-4a08-a805-3628048ebf1c.rna_seq.augmented_star_gene_counts.tsv
    ## 6 3ffa2cdf-4a54-4362-a6cf-b0567ef69551.rna_seq.augmented_star_gene_counts.tsv
    ##                           submitter_id           data_category            type
    ## 1 162ed5c3-add6-4485-8d99-37c02751bd36 Transcriptome Profiling gene_expression
    ## 2 95f3a1dc-5651-4ff8-965b-3d6820385cf0 Transcriptome Profiling gene_expression
    ## 3 4a83ecea-8c9c-4c2b-a753-f2479410dfef Transcriptome Profiling gene_expression
    ## 4 d182c7d9-0a71-4a53-9c0a-78ec0b5aad7e Transcriptome Profiling gene_expression
    ## 5 bef105f5-9f26-4122-bf53-762a2c5243e1 Transcriptome Profiling gene_expression
    ## 6 85fdc9b8-8383-429e-88b2-e3c170d651b4 Transcriptome Profiling gene_expression
    ##   file_size                 created_datetime                           md5sum
    ## 1   4220759 2021-12-13T22:37:49.340267-06:00 cae9cf46552fd8b701916e21c6663ad2
    ## 2   4242201 2021-12-13T22:05:24.743515-06:00 2edc4071ee61ac47fd0ff24e61e353d4
    ## 3   4240106 2021-12-13T22:33:53.775781-06:00 24f2aa8e80fb2cdae453f315f388cc09
    ## 4   4248363 2021-12-13T22:14:55.917966-06:00 2eefcb6236111244626b39ec34b8ed09
    ## 5   4249633 2021-12-13T21:51:47.586584-06:00 dd0fdbda13c78711977c48efacbd3026
    ## 6   4247939 2021-12-13T22:05:29.355511-06:00 1b927939b68699db0a1431d1da33007c
    ##                   updated_datetime                              file_id
    ## 1 2022-01-19T12:00:30.311996-06:00 a5c799bc-1332-4127-8fda-84da04084d25
    ## 2 2022-01-19T12:00:30.311996-06:00 1dd898ab-3b85-4932-b889-0eb0a494f2b6
    ## 3 2022-01-19T12:00:30.311996-06:00 7286a685-17ef-4ef5-9b7b-2696c6e7e1a6
    ## 4 2022-01-19T12:00:30.311996-06:00 de6f1503-33d7-4f86-b835-bdffba7ea4e3
    ## 5 2022-01-19T12:00:30.311996-06:00 4c46a463-f0e3-4e19-ad95-5dab0231d1c6
    ## 6 2022-01-19T12:00:30.311996-06:00 586d9c21-e58a-4171-9755-30c0e8b16293
    ##                        data_type    state experimental_strategy version
    ## 1 Gene Expression Quantification released               RNA-Seq       1
    ## 2 Gene Expression Quantification released               RNA-Seq       1
    ## 3 Gene Expression Quantification released               RNA-Seq       1
    ## 4 Gene Expression Quantification released               RNA-Seq       1
    ## 5 Gene Expression Quantification released               RNA-Seq       1
    ## 6 Gene Expression Quantification released               RNA-Seq       1
    ##   data_release   project                          analysis_id analysis_state
    ## 1  32.0 - 38.0 TCGA-BRCA 1d06d9a6-4145-4c38-a35d-c23c3433423a       released
    ## 2  32.0 - 38.0 TCGA-BRCA 0a4376f3-7132-4bf1-86bf-df060886982c       released
    ## 3  32.0 - 38.0 TCGA-BRCA d818e0d0-6685-4b7c-9575-d7f777f8943f       released
    ## 4  32.0 - 38.0 TCGA-BRCA 7603e12d-4e87-4225-b928-3c6e41effb28       released
    ## 5  32.0 - 38.0 TCGA-BRCA 339fee09-9a7e-4ddf-90e1-9144f22d1369       released
    ## 6  32.0 - 38.0 TCGA-BRCA c0a48ae9-c8b6-4a1d-8fd5-66c45399acbb       released
    ##                               analysis_submitter_id
    ## 1 de01516e-43f0-4f96-8ac6-ab543a314829_star__counts
    ## 2 e0e050cc-5e2a-41fc-8394-b4816a9dc8a2_star__counts
    ## 3 3572211c-0611-4d67-9e71-c5a8f704bd1f_star__counts
    ## 4 195d8bd4-352c-4407-9319-797d915c72b8_star__counts
    ## 5 f2ab4258-222d-4a08-a805-3628048ebf1c_star__counts
    ## 6 3ffa2cdf-4a54-4362-a6cf-b0567ef69551_star__counts
    ##                                                                                                                         analysis_workflow_link
    ## 1 https://github.com/NCI-GDC/gdc-rnaseq-cwl/blob/5d8c131bbff59fb0c969217fc1d44e6d1503cd1f/rnaseq-star-align/star2pass.rnaseq_harmonization.cwl
    ## 2 https://github.com/NCI-GDC/gdc-rnaseq-cwl/blob/5d8c131bbff59fb0c969217fc1d44e6d1503cd1f/rnaseq-star-align/star2pass.rnaseq_harmonization.cwl
    ## 3 https://github.com/NCI-GDC/gdc-rnaseq-cwl/blob/5d8c131bbff59fb0c969217fc1d44e6d1503cd1f/rnaseq-star-align/star2pass.rnaseq_harmonization.cwl
    ## 4 https://github.com/NCI-GDC/gdc-rnaseq-cwl/blob/5d8c131bbff59fb0c969217fc1d44e6d1503cd1f/rnaseq-star-align/star2pass.rnaseq_harmonization.cwl
    ## 5 https://github.com/NCI-GDC/gdc-rnaseq-cwl/blob/5d8c131bbff59fb0c969217fc1d44e6d1503cd1f/rnaseq-star-align/star2pass.rnaseq_harmonization.cwl
    ## 6 https://github.com/NCI-GDC/gdc-rnaseq-cwl/blob/5d8c131bbff59fb0c969217fc1d44e6d1503cd1f/rnaseq-star-align/star2pass.rnaseq_harmonization.cwl
    ##   analysis_workflow_type                analysis_workflow_version   sample_type
    ## 1          STAR - Counts 5d8c131bbff59fb0c969217fc1d44e6d1503cd1f Primary Tumor
    ## 2          STAR - Counts 5d8c131bbff59fb0c969217fc1d44e6d1503cd1f Primary Tumor
    ## 3          STAR - Counts 5d8c131bbff59fb0c969217fc1d44e6d1503cd1f Primary Tumor
    ## 4          STAR - Counts 5d8c131bbff59fb0c969217fc1d44e6d1503cd1f Primary Tumor
    ## 5          STAR - Counts 5d8c131bbff59fb0c969217fc1d44e6d1503cd1f Primary Tumor
    ## 6          STAR - Counts 5d8c131bbff59fb0c969217fc1d44e6d1503cd1f Primary Tumor
    ##   is_ffpe cases.submitter_id sample.submitter_id
    ## 1      NA       TCGA-A7-A26E    TCGA-A7-A26E-01B
    ## 2      NA       TCGA-A2-A0CU    TCGA-A2-A0CU-01A
    ## 3      NA       TCGA-PL-A8LV    TCGA-PL-A8LV-01A
    ## 4      NA       TCGA-BH-A0BC    TCGA-BH-A0BC-01A
    ## 5      NA       TCGA-AR-A1AX    TCGA-AR-A1AX-01A
    ## 6      NA       TCGA-AC-A2FO    TCGA-AC-A2FO-01A

``` r
normalMatched<-primary_query[[1]][[1]] %>% filter( sample_type == "Solid Tissue Normal")
```

``` r
normalMatched$cases[1]
```

    ## [1] "TCGA-BH-A1FU-11A-23R-A14D-07"

``` r
which(substr(normalMatched$cases,start = 1, stop = 12) %in% sample_data$BARCODE)
```

    ##  [1]  18  25  26  41  49  52  59  67  77  82  89  94  96 106 107 109

``` r
  ids<-normalMatched%>%slice(which(substr(normalMatched$cases,start = 1, stop = 12) %in% sample_data$BARCODE))%>%
  select(cases)%>%
  unlist(,use.names = F)
```

``` r
ret_q<-GDCquery(project = "TCGA-BRCA",data.category = "Transcriptome Profiling",
         data.type = "Gene Expression Quantification",
         workflow.type = "STAR - Counts",
         sample.type = c("Solid Tissue Normal"),
        barcode = ids
         )
```

    ## --------------------------------------

    ## o GDCquery: Searching in GDC database

    ## --------------------------------------

    ## Genome of reference: hg38

    ## --------------------------------------------

    ## oo Accessing GDC. This might take a while...

    ## --------------------------------------------

    ## ooo Project: TCGA-BRCA

    ## --------------------

    ## oo Filtering results

    ## --------------------

    ## ooo By data.type

    ## ooo By workflow.type

    ## ooo By barcode

    ## ooo By sample.type

    ## ----------------

    ## oo Checking data

    ## ----------------

    ## ooo Checking if there are duplicated cases

    ## ooo Checking if there are results for the query

    ## -------------------

    ## o Preparing output

    ## -------------------

``` r
GDCdownload(ret_q,files.per.chunk = 10)
```

    ## Downloading data for project TCGA-BRCA

    ## Of the 16 files for download 16 already exist.

    ## All samples have been already downloaded

``` r
dir1<-"./GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification"

files_list<-list.files(dir1, full.names = T,recursive = T)
```

``` r
file_names<-str_replace(dir(dir1,full.names = F, recursive = T), ".*\\/(.*)\\.rna_seq.*", "\\1")

samples_list<-lapply(files_list, function(x){
  t1<-read.table(x,sep = "\t", header = T)%>%
    filter(gene_type == "protein_coding")%>%
    select(c(4))
  return(t1)
})
```

``` r
a<-do.call(base::cbind,samples_list)
```

``` r
colnames(a)<-unlist(file_names)
```

``` r
normal_df<-data.frame(barcode = substr(ids, start = 1, stop = 12), file_name_normal = file_names)
```

``` r
head(normal_df)
```

    ##        barcode                     file_name_normal
    ## 1 TCGA-BH-A0B3 d4f91697-1c39-4398-bf2f-85217a22ddff
    ## 2 TCGA-E2-A1LH 678aa892-0631-47cc-b19e-cdcef42cecb3
    ## 3 TCGA-BH-A18Q 930d7c02-c230-4d19-9247-74fb0e02d524
    ## 4 TCGA-E9-A1RH 308e15c9-16b0-4f24-b9dc-e9dc91038579
    ## 5 TCGA-BH-A0BW d23ff3bc-5524-4dca-b0f4-a1561a30566c
    ## 6 TCGA-BH-A0E0 0a950719-ca7a-4dbf-8755-262cabf3d4b7

``` r
case_df<-sample_data%>%select(BARCODE, file_name,TNBCtype,TNBCtype_4,PAM50,PAM50lite)%>%
mutate(file_name = str_replace(file_name ,"(.*)\\.rna_seq.*", "\\1" ) )%>%
 slice(which(BARCODE %in% normal_df$barcode))
```

``` r
head(case_df)
```

    ##        BARCODE                            file_name TNBCtype TNBCtype_4  PAM50
    ## 1 TCGA-BH-A0B3 41bdc474-8e6e-4f39-98ec-6d73ee41837c      UNC        LAR   LumA
    ## 2 TCGA-BH-A0E0 2b5864c9-7577-4e93-8d22-1f45a98267f1      LAR        LAR   LumA
    ## 3 TCGA-BH-A1F6 1a2bc49d-edc6-4b84-88f0-73b6ed250d1b      MSL        LAR   LumA
    ## 4 TCGA-BH-A1FC 31c69703-e318-4e31-aae8-c6cfabe61bc3      MSL        LAR   LumA
    ## 5 TCGA-GI-A2C9 b161c26d-64b1-4e03-94da-4f9b9de25e88      MSL          M Normal
    ## 6 TCGA-E2-A1L7 2056c543-2f73-4a2f-995e-3af2dda7b514      MSL        LAR Normal
    ##   PAM50lite
    ## 1 Non-basal
    ## 2 Non-basal
    ## 3 Non-basal
    ## 4 Non-basal
    ## 5 Non-basal
    ## 6 Non-basal

``` r
colnames(case_df)[1]<-"barcode"
```

``` r
matchedCaseControl<-inner_join(case_df, normal_df,"barcode")
```

``` r
head(matchedCaseControl)
```

    ##        barcode                            file_name TNBCtype TNBCtype_4  PAM50
    ## 1 TCGA-BH-A0B3 41bdc474-8e6e-4f39-98ec-6d73ee41837c      UNC        LAR   LumA
    ## 2 TCGA-BH-A0E0 2b5864c9-7577-4e93-8d22-1f45a98267f1      LAR        LAR   LumA
    ## 3 TCGA-BH-A1F6 1a2bc49d-edc6-4b84-88f0-73b6ed250d1b      MSL        LAR   LumA
    ## 4 TCGA-BH-A1FC 31c69703-e318-4e31-aae8-c6cfabe61bc3      MSL        LAR   LumA
    ## 5 TCGA-GI-A2C9 b161c26d-64b1-4e03-94da-4f9b9de25e88      MSL          M Normal
    ## 6 TCGA-E2-A1L7 2056c543-2f73-4a2f-995e-3af2dda7b514      MSL        LAR Normal
    ##   PAM50lite                     file_name_normal
    ## 1 Non-basal d4f91697-1c39-4398-bf2f-85217a22ddff
    ## 2 Non-basal 0a950719-ca7a-4dbf-8755-262cabf3d4b7
    ## 3 Non-basal 7be7033e-db49-45cd-8ea5-827266707b1e
    ## 4 Non-basal a4cf670f-1627-42cf-9917-fc2c89708353
    ## 5 Non-basal 898f9214-8dab-4b92-9233-70c06c60078c
    ## 6 Non-basal da752364-731c-4e77-8d6e-9c91a5b7f68f

``` r
colnames(matchedCaseControl)[2]<-"file_name_tnbc"
```

``` r
ex1<-read.table(files_list[1],sep = "\t", header = T)%>%
  filter(gene_type == "protein_coding")%>%
  select(c(1,2,4))
```

``` r
rownames(a)<-ex1$gene_id
```

``` r
a<-a[!grepl("_PAR_Y",ex1$gene_id),]
```

``` r
rownames(a)<-str_replace(rownames(a), "^(ENSG\\d+)[:punct:]\\d+","\\1" )
```

``` r
library(readr)
tnbc_counts<-read_csv("./samples_from_tcga/count_file_names_final.csv")
```

    ## New names:
    ## Rows: 19944 Columns: 184
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): ...1 dbl (183): e0527b4b-de23-4b9f-8871-13f64f6976cd,
    ## 8a27980e-a506-4cb3-91f2-3f0...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
tnbc_counts<-tnbc_counts%>%as.data.frame()%>%
  set_rownames(.$...1)%>%dplyr::select(-c(1))
```

``` r
mathcedtnbc<-tnbc_counts%>%
  select(which(colnames(tnbc_counts)%in%matchedCaseControl$file_name_tnbc))
```

``` r
countsMatched<-cbind(mathcedtnbc,a)
```

``` r
head(countsMatched)
```

    ##                 41bdc474-8e6e-4f39-98ec-6d73ee41837c
    ## ENSG00000000003                                 2480
    ## ENSG00000000005                                   11
    ## ENSG00000000419                                 1845
    ## ENSG00000000457                                  949
    ## ENSG00000000460                                 1024
    ## ENSG00000000938                                  803
    ##                 2b5864c9-7577-4e93-8d22-1f45a98267f1
    ## ENSG00000000003                                 2836
    ## ENSG00000000005                                   33
    ## ENSG00000000419                                 1886
    ## ENSG00000000457                                 1473
    ## ENSG00000000460                                  937
    ## ENSG00000000938                                  942
    ##                 1a2bc49d-edc6-4b84-88f0-73b6ed250d1b
    ## ENSG00000000003                                 3954
    ## ENSG00000000005                                    4
    ## ENSG00000000419                                 2591
    ## ENSG00000000457                                  565
    ## ENSG00000000460                                  413
    ## ENSG00000000938                                  547
    ##                 31c69703-e318-4e31-aae8-c6cfabe61bc3
    ## ENSG00000000003                                 4730
    ## ENSG00000000005                                    2
    ## ENSG00000000419                                 2284
    ## ENSG00000000457                                 1549
    ## ENSG00000000460                                 1510
    ## ENSG00000000938                                  760
    ##                 b161c26d-64b1-4e03-94da-4f9b9de25e88
    ## ENSG00000000003                                 3926
    ## ENSG00000000005                                    2
    ## ENSG00000000419                                 2696
    ## ENSG00000000457                                 1271
    ## ENSG00000000460                                  585
    ## ENSG00000000938                                  418
    ##                 2056c543-2f73-4a2f-995e-3af2dda7b514
    ## ENSG00000000003                                 1689
    ## ENSG00000000005                                   16
    ## ENSG00000000419                                 1810
    ## ENSG00000000457                                 1098
    ## ENSG00000000460                                  715
    ## ENSG00000000938                                  624
    ##                 d0ee5ff7-a49a-4633-93a6-40c9e29fb0b7
    ## ENSG00000000003                                  324
    ## ENSG00000000005                                    0
    ## ENSG00000000419                                 6231
    ## ENSG00000000457                                  661
    ## ENSG00000000460                                  695
    ## ENSG00000000938                                  333
    ##                 58b399a0-3070-44aa-9a40-f0e2c0fea0fc
    ## ENSG00000000003                                 8431
    ## ENSG00000000005                                   54
    ## ENSG00000000419                                 2399
    ## ENSG00000000457                                 1292
    ## ENSG00000000460                                 1845
    ## ENSG00000000938                                  284
    ##                 f48956a3-c821-429f-b1ad-a27133d739bd
    ## ENSG00000000003                                 3600
    ## ENSG00000000005                                   77
    ## ENSG00000000419                                 2078
    ## ENSG00000000457                                 2576
    ## ENSG00000000460                                  948
    ## ENSG00000000938                                  349
    ##                 9440b088-b04e-4417-9543-5001b0a2b217
    ## ENSG00000000003                                 1091
    ## ENSG00000000005                                   24
    ## ENSG00000000419                                 3559
    ## ENSG00000000457                                 1425
    ## ENSG00000000460                                 1967
    ## ENSG00000000938                                  847
    ##                 5e4626c1-dcf9-4781-a2e7-f6afca7e2ee7
    ## ENSG00000000003                                 2690
    ## ENSG00000000005                                   27
    ## ENSG00000000419                                 2435
    ## ENSG00000000457                                 2226
    ## ENSG00000000460                                  945
    ## ENSG00000000938                                  571
    ##                 28b12d43-f858-458b-9305-6581a9498e59
    ## ENSG00000000003                                  680
    ## ENSG00000000005                                   24
    ## ENSG00000000419                                 1184
    ## ENSG00000000457                                  559
    ## ENSG00000000460                                  254
    ## ENSG00000000938                                  881
    ##                 7880b3ba-9f39-4744-be72-b87dbd8eb425
    ## ENSG00000000003                                 3786
    ## ENSG00000000005                                    1
    ## ENSG00000000419                                 1639
    ## ENSG00000000457                                  787
    ## ENSG00000000460                                  568
    ## ENSG00000000938                                  893
    ##                 dd1d9d28-3929-4361-a6ff-0ac781fb0ac3
    ## ENSG00000000003                                 8037
    ## ENSG00000000005                                    1
    ## ENSG00000000419                                 2814
    ## ENSG00000000457                                 1920
    ## ENSG00000000460                                  795
    ## ENSG00000000938                                  163
    ##                 8bf23a6b-8ca6-4fcf-b4e9-a14c3a875526
    ## ENSG00000000003                                 6693
    ## ENSG00000000005                                    6
    ## ENSG00000000419                                 2246
    ## ENSG00000000457                                 3356
    ## ENSG00000000460                                  401
    ## ENSG00000000938                                   93
    ##                 aeee03c3-20b3-4f96-b686-52111ce99a63
    ## ENSG00000000003                                 3211
    ## ENSG00000000005                                    1
    ## ENSG00000000419                                 2144
    ## ENSG00000000457                                 1707
    ## ENSG00000000460                                  522
    ## ENSG00000000938                                  571
    ##                 d4f91697-1c39-4398-bf2f-85217a22ddff
    ## ENSG00000000003                                 3116
    ## ENSG00000000005                                  817
    ## ENSG00000000419                                 1408
    ## ENSG00000000457                                  954
    ## ENSG00000000460                                  125
    ## ENSG00000000938                                  461
    ##                 678aa892-0631-47cc-b19e-cdcef42cecb3
    ## ENSG00000000003                                 4468
    ## ENSG00000000005                                  109
    ## ENSG00000000419                                 1291
    ## ENSG00000000457                                 1392
    ## ENSG00000000460                                  267
    ## ENSG00000000938                                  406
    ##                 930d7c02-c230-4d19-9247-74fb0e02d524
    ## ENSG00000000003                                 3745
    ## ENSG00000000005                                  254
    ## ENSG00000000419                                 1582
    ## ENSG00000000457                                 1003
    ## ENSG00000000460                                  177
    ## ENSG00000000938                                  184
    ##                 308e15c9-16b0-4f24-b9dc-e9dc91038579
    ## ENSG00000000003                                 5383
    ## ENSG00000000005                                  319
    ## ENSG00000000419                                 2465
    ## ENSG00000000457                                 1191
    ## ENSG00000000460                                  272
    ## ENSG00000000938                                  862
    ##                 d23ff3bc-5524-4dca-b0f4-a1561a30566c
    ## ENSG00000000003                                 3051
    ## ENSG00000000005                                 1083
    ## ENSG00000000419                                 1359
    ## ENSG00000000457                                  755
    ## ENSG00000000460                                  182
    ## ENSG00000000938                                  846
    ##                 0a950719-ca7a-4dbf-8755-262cabf3d4b7
    ## ENSG00000000003                                 7631
    ## ENSG00000000005                                  191
    ## ENSG00000000419                                 2150
    ## ENSG00000000457                                 2237
    ## ENSG00000000460                                  568
    ## ENSG00000000938                                  572
    ##                 d83bb114-c175-4fd5-aac0-5dc04b9218b0
    ## ENSG00000000003                                 4816
    ## ENSG00000000005                                  308
    ## ENSG00000000419                                 2279
    ## ENSG00000000457                                 1086
    ## ENSG00000000460                                  294
    ## ENSG00000000938                                  967
    ##                 d920c232-82ae-48f9-98d1-b94bf4a45d8d
    ## ENSG00000000003                                 6764
    ## ENSG00000000005                                  462
    ## ENSG00000000419                                 1916
    ## ENSG00000000457                                 1888
    ## ENSG00000000460                                  377
    ## ENSG00000000938                                  567
    ##                 7be7033e-db49-45cd-8ea5-827266707b1e
    ## ENSG00000000003                                 2423
    ## ENSG00000000005                                  875
    ## ENSG00000000419                                 1244
    ## ENSG00000000457                                  554
    ## ENSG00000000460                                  163
    ## ENSG00000000938                                 4942
    ##                 269c35f0-a4f7-4e30-a69f-f1f3b7b5dace
    ## ENSG00000000003                                 3573
    ## ENSG00000000005                                  770
    ## ENSG00000000419                                 1607
    ## ENSG00000000457                                 1040
    ## ENSG00000000460                                  215
    ## ENSG00000000938                                  573
    ##                 a4cf670f-1627-42cf-9917-fc2c89708353
    ## ENSG00000000003                                 3994
    ## ENSG00000000005                                 1508
    ## ENSG00000000419                                 1963
    ## ENSG00000000457                                  903
    ## ENSG00000000460                                  246
    ## ENSG00000000938                                  412
    ##                 484b85c7-1102-4d34-89b8-8ce458b55d2b
    ## ENSG00000000003                                 1816
    ## ENSG00000000005                                  682
    ## ENSG00000000419                                  875
    ## ENSG00000000457                                  489
    ## ENSG00000000460                                  149
    ## ENSG00000000938                                  502
    ##                 898f9214-8dab-4b92-9233-70c06c60078c
    ## ENSG00000000003                                 7064
    ## ENSG00000000005                                 1899
    ## ENSG00000000419                                 2680
    ## ENSG00000000457                                 1626
    ## ENSG00000000460                                  527
    ## ENSG00000000938                                 1235
    ##                 da752364-731c-4e77-8d6e-9c91a5b7f68f
    ## ENSG00000000003                                 4209
    ## ENSG00000000005                                   71
    ## ENSG00000000419                                 1611
    ## ENSG00000000457                                 1217
    ## ENSG00000000460                                  346
    ## ENSG00000000938                                  787
    ##                 078f6959-adcf-44a5-91d2-2581ba3c8e62
    ## ENSG00000000003                                10613
    ## ENSG00000000005                                 1465
    ## ENSG00000000419                                 4030
    ## ENSG00000000457                                 3392
    ## ENSG00000000460                                 1082
    ## ENSG00000000938                                  606
    ##                 741fb6d6-364f-47a0-a262-f18235c48fc4
    ## ENSG00000000003                                 2690
    ## ENSG00000000005                                   52
    ## ENSG00000000419                                 1529
    ## ENSG00000000457                                 1271
    ## ENSG00000000460                                  388
    ## ENSG00000000938                                 1638
