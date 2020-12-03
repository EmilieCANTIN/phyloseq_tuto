Controle\_continu\_1\_analyse\_des\_donnees
================

# Méthodes

## De lectures brutes aux tableaux

Ce premier code permet d’importer les données d’une étude Illumina
MiSeq, à partir d’un ensemble de fichiers fastq.Ici, on définit une
variable chemin miseq\_path, afin de pouvoir accéder à ces données.

``` r
library("dada2")
```

    ## Loading required package: Rcpp

``` r
miseq_path <- "~/MiSeq_SOP" # MODIFIER le répertoire contenant les fichiers fastq après la décompression.
list.files(miseq_path)
```

    ##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
    ##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
    ##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
    ##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
    ##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
    ## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
    ## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
    ## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
    ## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
    ## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
    ## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
    ## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
    ## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
    ## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
    ## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
    ## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
    ## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
    ## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
    ## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
    ## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
    ## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
    ## [43] "mouse.dpw.metadata"            "mouse.time.design"            
    ## [45] "stability.batch"               "stability.files"

## Filtrer les données

On filtre ensuite les séquences de faible qualité, puis on les enlève.
On demande ici d’afficher les “moins bons”.

``` r
# Le tri permet de s'assurer que les lectures en avant et en arrière sont dans le même ordre
fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
# Extraire les noms des échantillons, en supposant que les noms de fichiers ont un format : SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Préciser le chemin complet vers les fnFs et fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
fnFs[1:3]
```

    ## [1] "~/MiSeq_SOP/F3D0_S188_L001_R1_001.fastq"  
    ## [2] "~/MiSeq_SOP/F3D1_S189_L001_R1_001.fastq"  
    ## [3] "~/MiSeq_SOP/F3D141_S207_L001_R1_001.fastq"

``` r
fnRs[1:3]
```

    ## [1] "~/MiSeq_SOP/F3D0_S188_L001_R2_001.fastq"  
    ## [2] "~/MiSeq_SOP/F3D1_S189_L001_R2_001.fastq"  
    ## [3] "~/MiSeq_SOP/F3D141_S207_L001_R2_001.fastq"

On sait que plus on se rapproche de la fin des séquençages, moins bonne
sera leur qualité. En effet, on remarque que pour les lectures avant
(deux premiers graphes), le score de qualité moyen ne descend jamais en
dessous de 25. Au contraire, les graphes incarnant la fin des lectures
montrent un score de qualité plus bas (\~20). ce type de chiffre
représente la probabilité que ce ne soit pas le bon nucléotide
d’appelé. De ce fait, avec un Q20 en fin de séquences, il y a une
chance sur 100 que ce soit le cas.

``` r
library("dada2")
library("ggplot2")
plotQualityProfile(fnFs[1:2])
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
plotQualityProfile(fnRs[1:2])
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
On voit bien d’après ces graphiques ci-dessus que les scores de qualités
baissent vers la position 240 pour les premières lectures, et plutôt
vers la position 160 pour les lectures arrières. En prenant ces
informations en compte, on va pouvoir dans un premier temps créer des
variables pour les fichiers filtrés, puis appliquer la fonction
filterAndTrim.

``` r
filt_path <- file.path(miseq_path, "filtered") # Placez les fichiers filtrés dans le sous-répertoire filtered/
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))
```

### Filtrez les lectures en amont et en aval

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
head(out)
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

Cette fonction se base sur des fichiers contenant les lectures coupées
ayant passées les filtres.

## Variantes de séquences d’inférence

### Dereplication

Ce type de fonction diminue sensiblement le temps de calcul des codes à
suivre, en supprimant les comparaisons redondantes.Les résultats
ressortants marquent le nombre de lectures à séquence unique, pour
chaque fichier.

``` r
derepFs <- derepFastq(filtFs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D0_F_filt.fastq.gz

    ## Encountered 1979 unique sequences from 7113 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D1_F_filt.fastq.gz

    ## Encountered 1639 unique sequences from 5299 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D141_F_filt.fastq.gz

    ## Encountered 1477 unique sequences from 5463 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D142_F_filt.fastq.gz

    ## Encountered 904 unique sequences from 2914 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D143_F_filt.fastq.gz

    ## Encountered 939 unique sequences from 2941 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D144_F_filt.fastq.gz

    ## Encountered 1267 unique sequences from 4312 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D145_F_filt.fastq.gz

    ## Encountered 1756 unique sequences from 6741 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D146_F_filt.fastq.gz

    ## Encountered 1438 unique sequences from 4560 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D147_F_filt.fastq.gz

    ## Encountered 3590 unique sequences from 15637 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D148_F_filt.fastq.gz

    ## Encountered 2762 unique sequences from 11413 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D149_F_filt.fastq.gz

    ## Encountered 3021 unique sequences from 12017 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D150_F_filt.fastq.gz

    ## Encountered 1566 unique sequences from 5032 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D2_F_filt.fastq.gz

    ## Encountered 3707 unique sequences from 18075 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D3_F_filt.fastq.gz

    ## Encountered 1479 unique sequences from 6250 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D5_F_filt.fastq.gz

    ## Encountered 1195 unique sequences from 4052 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D6_F_filt.fastq.gz

    ## Encountered 1832 unique sequences from 7369 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D7_F_filt.fastq.gz

    ## Encountered 1183 unique sequences from 4765 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D8_F_filt.fastq.gz

    ## Encountered 1382 unique sequences from 4871 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D9_F_filt.fastq.gz

    ## Encountered 1709 unique sequences from 6504 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/Mock_F_filt.fastq.gz

    ## Encountered 897 unique sequences from 4314 total sequences read.

``` r
derepRs <- derepFastq(filtRs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D0_R_filt.fastq.gz

    ## Encountered 1660 unique sequences from 7113 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D1_R_filt.fastq.gz

    ## Encountered 1349 unique sequences from 5299 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D141_R_filt.fastq.gz

    ## Encountered 1335 unique sequences from 5463 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D142_R_filt.fastq.gz

    ## Encountered 853 unique sequences from 2914 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D143_R_filt.fastq.gz

    ## Encountered 880 unique sequences from 2941 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D144_R_filt.fastq.gz

    ## Encountered 1286 unique sequences from 4312 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D145_R_filt.fastq.gz

    ## Encountered 1803 unique sequences from 6741 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D146_R_filt.fastq.gz

    ## Encountered 1265 unique sequences from 4560 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D147_R_filt.fastq.gz

    ## Encountered 3414 unique sequences from 15637 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D148_R_filt.fastq.gz

    ## Encountered 2522 unique sequences from 11413 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D149_R_filt.fastq.gz

    ## Encountered 2771 unique sequences from 12017 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D150_R_filt.fastq.gz

    ## Encountered 1415 unique sequences from 5032 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D2_R_filt.fastq.gz

    ## Encountered 3290 unique sequences from 18075 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D3_R_filt.fastq.gz

    ## Encountered 1390 unique sequences from 6250 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D5_R_filt.fastq.gz

    ## Encountered 1134 unique sequences from 4052 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D6_R_filt.fastq.gz

    ## Encountered 1635 unique sequences from 7369 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D7_R_filt.fastq.gz

    ## Encountered 1084 unique sequences from 4765 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D8_R_filt.fastq.gz

    ## Encountered 1161 unique sequences from 4871 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/F3D9_R_filt.fastq.gz

    ## Encountered 1502 unique sequences from 6504 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/MiSeq_SOP/filtered/Mock_R_filt.fastq.gz

    ## Encountered 732 unique sequences from 4314 total sequences read.

``` r
# Nommer les objets de la classe derep par les noms des échantillons
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames
```

Ci-dessous, la fonction learnErrors permet d’estimer les taux d’erreurs
à partir d’un grand ensemble de données. Ainsi, les résultats ci-après
expriment le nombre de bases qui sera finalement utilisé, par rapport au
premier ensemble.

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
library("dada2")
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
library("dada2")
plotErrors(errR, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->
Les figures ci-dessus représentent les estimations des taux d’erreurs.
La ligne rouge incarne la tendance générale du graphique. Ensuite, les
points noirs reflètent le taux d’erreurs observées, et la ligne noire le
taux d’erreurs ajustées. On peut donc observer ci-dessus la fréquence du
taux d’erreur en fonction du score de qualité. Aucune différence
significative ne peut être relevée entre errR et errF. En effet, on
observe la même tendance : moins il y a d’erreurs, plus le score de
qualité augmente, ce qui est en accord avec les résultats attendus.

``` r
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

``` r
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

Les résultats ci-dessus signifient que 128 séquences ont été
spécialement extraites et définies comme des variantes réelles. Elles
ont été déterminées à partir d’un ensemble de 1979 séquences uniques.

## Construire un tableau de séquences et éliminer les chimères

``` r
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
table(nchar(getSequences(seqtabAll)))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  85 186   5   2

Ici,les chimères n’ont pas encore été enlevées. Le code ci-dessous y
remédie. En effet, il supprime les séquences reproduites en comparant
chaque séquence aux autres.

``` r
seqtabNoC <- removeBimeraDenovo(seqtabAll)
dim(seqtabNoC)
```

    ## [1]  19 218

On peut donc dire que les chimères représentent environ 19% des
variantes de séquences.

## Attribuer une taxonomie

On va chercher ici à classer toutes les variantes des séquences
étudiées. Cela est rendu possible car l’ARN 16S est un marqueur
extrêmement bien classé. Tout d’abord, il nous faut importer les
données Silva, qui nous serviront à réaliser l’arbre taxonomique. Vous
pourrez retrouver les codes nécessaires dans 01\_data\_import.

``` r
taxa <- assignTaxonomy(seqtabNoC, "~/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
```

``` r
taxa.print <- taxa # Suppression des noms de séquences pour l'affichage uniquement
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum         Class         Order           Family          
    ## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      Genus        
    ## [1,] NA           
    ## [2,] NA           
    ## [3,] NA           
    ## [4,] NA           
    ## [5,] "Bacteroides"
    ## [6,] NA

On peut ici remarquer que les plus représentés sont les bactéroides.
Cela est peu surprenant car ils représentent une partie importante de la
flore intestinale.

# PhyloSeq

Avant toute chose, il est nécessaire d’installer quelques packages pour
commencer. Vous pourrez les retrouver dans 00\_package\_installation.
Par exemple, nous aurons besoin de Decipher (ou Bioconductor) - permet
de déchiffrer et gérer des données - et de Phangorn (ou Cran). Ensuite,
les lignes suivantes permettent de créer des varibles nécessaires à la
construction d’un arbre phylogénétique.

``` r
samples.out <- rownames(seqtabNoC)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out # 
```

## Combiner des données dans un objet phyloseq

Dès lors, les données utilisées jusque là peuvent être combinées. C’est
à dire que les séquences, la taxonomie, et notre arbre formeront un
même objet. Ici, nous allons créer une variable ps.

``` r
library(phangorn)
```

    ## Loading required package: ape

``` r
library(DECIPHER)
```

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:ape':
    ## 
    ##     complement

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: RSQLite

``` r
seqs <- getSequences(seqtabNoC)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
```

    ## negative edges length changed to 0!

``` r
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
        rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
```

``` r
library(phyloseq)
```

    ## 
    ## Attaching package: 'phyloseq'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     distance

``` r
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa),phy_tree(fitGTR$tree))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Retirer l'échantillon fictif
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 218 taxa and 19 samples ]
    ## sample_data() Sample Data:       [ 19 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 218 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 218 tips and 216 internal nodes ]

## Chargement des données

``` r
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 218 taxa and 19 samples ]
    ## sample_data() Sample Data:       [ 19 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 218 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 218 tips and 216 internal nodes ]
    ## refseq()      DNAStringSet:      [ 218 reference sequences ]

## Filtrage

Cette étape est notamment réalisée pour éviter de perdre du temps à
analyser un grand nombre de taxons observés que certaines fois (taxons
rares).

### Filtrage taxonomique

``` r
ps_connect <-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/ps.rds") 
ps = readRDS(ps_connect)
```

``` r
# Afficher les rangs disponibles dans l'ensemble de données
rank_names(ps)
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"

``` r
# Créer un tableau, avec nombre de caractéristiques pour chaque phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)
```

    ## 
    ##              Actinobacteria               Bacteroidetes 
    ##                          13                          23 
    ## Candidatus_Saccharibacteria   Cyanobacteria/Chloroplast 
    ##                           1                           4 
    ##         Deinococcus-Thermus                  Firmicutes 
    ##                           1                         327 
    ##                Fusobacteria              Proteobacteria 
    ##                           1                          11 
    ##                 Tenericutes             Verrucomicrobia 
    ##                           1                           1 
    ##                        <NA> 
    ##                           6

On peut ici voir que une seule caractéristique a été relevée pour
Deinococcus, Candidatus, Fusobacteria, Tenericutes et Verrucomicrobia.
Il est donc possible que chacune d’entre elles soit filtrée. Cela va
être vérifié avec les codes ci-après.

``` r
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
```

Ici, nous nous sommes assurés que les séquences incertaines (NA) soient
soustraites. Ensuite, la variable prevdf incarne la prévalence de
l’ensemble de nos données. Elle permettra ensuite de calculer la
prévalence totale puis moyenne de chaque phylum observé.

``` r
# Calculer la prévalence de chaque caractéristique, puis la stocker sous forme de data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Ajoutez la taxonomie et le nombre total de lectures à ces data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
```

``` r
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
```

    ##                         Phylum         1     2
    ## 1               Actinobacteria 120.15385  1562
    ## 2                Bacteroidetes 265.52174  6107
    ## 3  Candidatus_Saccharibacteria 280.00000   280
    ## 4    Cyanobacteria/Chloroplast  64.25000   257
    ## 5          Deinococcus-Thermus  52.00000    52
    ## 6                   Firmicutes 179.24771 58614
    ## 7                 Fusobacteria   2.00000     2
    ## 8               Proteobacteria  59.09091   650
    ## 9                  Tenericutes 234.00000   234
    ## 10             Verrucomicrobia 104.00000   104

On peut donc ici observer la prévalence totale (colonne 2) et la
prévalence moyenne (colonne 1). Les deux phylum que l’on va donc
considérer comme rare sont Fusobacteria et Deinococcus. Ainsi, il va
être nécessaire de les filtrer.

``` r
# Définir les phyla à filtrer
filterPhyla = c("Fusobacteria", "Deinococcus-Thermus")
# Filtrer les entrées avec un phylum non identifié
ps1 = subset_taxa(ps, !Phylum %in% filterPhyla)
ps1
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 381 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 381 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 381 tips and 379 internal nodes ]

Les chiffres représentés ci-dessus sont donc exempts des fusobactéries
ainsi que de Deinococcus.

### Filtrage de la prévalence

``` r
# Sous-ensemble du reste du phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Inclure une estimation pour le paramètre
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->
On peut ici observer la relation entre la prévalence et le nombre total
de lecture pour chaque phyla. Les firmicutes sont très présents quand il
y a très peu de Candidatus, Tenericutes ou encore Verrucomicrobia. Sur
chacun des graphiques on peut relever une certaine tendance même si on
remarque tout de même la présence de valeurs aberrantes (qu’il pourrait
être judicieux d’enlever pour des futures études). Ainsi, des
graphiques de ce type peuvent être particulièrement utiles pour
déterminer des caractéristiques de filtrage. On peut d’aileurs
également y remarquer une ligne en pointillé incarant le seuil de
prévalence utilisé ici.

``` r
# Définir le seuil de prévalence à 5 % du total des échantillons
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold
```

    ## [1] 18

``` r
# Exécuter le filtre de prévalence, en utilisant la fonction `prune_taxa()`
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)
```

## Agglomération taxonomique

``` r
# Combien de genres seraient présents après filtrage ?
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
```

    ## [1] 49

``` r
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)
```

``` r
h1 = 0.4
ps4 = tip_glom(ps2, h = h1)
```

Au lieu de passer par le rang taxonomique, il est préféré ici de définir
une hauteur d’arbre. Elle correspondra autrement dit à la distance
phylogénétique entre les différentes caractéristiques.

``` r
multiPlotTitleTextSize = 15
p2tree = plot_tree(ps2, method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(ps3, method = "treeonly",
                   ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p4tree = plot_tree(ps4, method = "treeonly",
                   ladderize = "left", title = "By Height") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
```

Le bloc de codes ci-dessus permet de comparer trois informations
différentes : données non filtrées avec arbre après agglomération
taxonomique et phylogénétique.

``` r
# Regroupement des 'plots'
library("gridExtra")
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
grid.arrange(nrow = 1, p2tree, p3tree, p4tree)
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->
Les arbres ci-dessus représentent donc différents types d’agglomération.
Dans un premier temps, on voit l’arbre original à gauche. Ensuite, celui
du milieu a été construit par agglomération par genre. Enfin, l’arbre
observé à droite incarne l’agglomération phylogénétique avec la distance
fixe définie ligne 300 (h1=0.4).

## Transformation de la valeur de l’abondance

``` r
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Sous-ensemble arbitraire, basé sur le Phylum, pour le tracé
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "sex",y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
```

``` r
# Transformer en abondance relative puis enregistrer comme nouvel objet
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})
```

``` r
plotBefore = plot_abundance(ps3,"")
plotAfter = plot_abundance(ps3ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 1,  plotBefore, plotAfter)
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->
Ces deux lots de graphiques permettent de comparer les abondances
initiales (à gauche) avec les données transformées (à droite). On peut
notamment relever la différence d’échelle entre les deux. Mise à part
cette différence (tout de même importante), les allures des “nuages de
points” restent sensiblement les mêmes. On peut notamment relever que la
distribution de Lactobacillus paraît être bimodale. Les graphiques
ci-après cherchent donc à le vérifier.

## Sous-ensemble par taxonomie

``` r
psOrd = subset_taxa(ps3ra, Order == "Lactobacillales")
plot_abundance(psOrd, Facet = "Genus", Color = NULL)
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->
Un sous-ensemble chez Lactobacillus a donc été crée. Ces graphiques
semblent représentés un tracé plus précis que les précédents. En effet,
on peut notamment y remarquer que la distribution des lactobacilles
semble être ici plutôt monomodale. On peut donc logiquement supposer que
les deux sous-ensembles représentés ici ont été mélangés. Ensuite, pour
réaliser des analyses plus poussées, nous aurons besoin d’un ensemble de
packages que vous pourrez retrouver dans cran (installés dans
00\_packages\_installation).

### Prétraitement

``` r
qplot(sample_data(ps)$age, geom = "histogram",binwidth=20) + xlab("age")
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->
Ces histogrammes représentent les différents compatages en fonction de
l’âge des souris. On peut notamment noter trois groupes différents qui
se distinguent. Plus les souris sont jeunes, plus le nombre de comptage
est élevé. Il pourrait donc être judicieux de relier l’âge des souris
aux pics de comptage.

``` r
qplot(log10(rowSums(otu_table(ps))),binwidth=0.2) +
  xlab("Logged counts-per-sample")
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->
Le graphique ci-dessus représente ainsi les comptages en fonction des
profondeurs de lectures log transformées. On peut voir, notamment par
l’allure de l’histogramme, que la transformation par logarithme n’est
pas suffisante pour normaliser en vue d’analyses. En effet, normalement,
la normalisation permet de rendre des ensembles comparables, or ici elle
n’est pas suffisante.

``` r
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship=gsub(" ","",sample_data(ps)$family_relationship)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGACTGCAAGTCAGATGTGAAACCCATGGGCTCAACCCATGGCCTGCATTTGAAACTGTAGTTCTTGAGTGATGGAGAGGCAGGCGGAATTCCGTGTGTAGCGGTGAAATGCGTAGATATACGGAGGAACACCAGTGGCGAAGGCGGCCTGCTGGACATTAACTGACGCTGAGGCGCGAAAGCGTGGGGAG
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->
Ce graphique représente une analyse par PCoA. On peut voir que les
souris on été classées par catégories d’âges. Par ailleurs, cette
analyse révèle notamment certaines erreurs parraissant aberrantes.

``` r
rel_abund <- t(apply(otu_table(ps), 1, function(x) x / sum(x)))
qplot(rel_abund[, 12], geom = "histogram",binwidth=0.05) +
  xlab("Relative abundance")
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->
Ici, les échantillons aberrants sont dominés par une seule VSA (=
variant de séquence d’amplicon). Etant un variant, on pourrait donc
s’attendre à ce que son abondance soit plutôt faible. Or ici,
l’abondance relative de cette catégorie est de largement plus de 90%,
même si dans l’ensemble sa diversité est très faible. En effet, il
s’agit bien là d’abondance relative et non absolue.

# Projections d’ordinations différentes

Désormais, nous allons enlever les valeurs aberrantes afin de pouvoir
mieux exploiter les données, notamment à travers des analyses par
ordination.

``` r
outliers <- c("F5D165", "F6D165", "M3D175", "M4D175", "M5D175", "M6D175")
ps <- prune_samples(!(sample_names(ps) %in% outliers), ps)
```

``` r
which(!rowSums(otu_table(ps)) > 1000)
```

    ## F5D145 M1D149   M1D9 M2D125  M2D19 M3D148 M3D149   M3D3   M3D5   M3D8 
    ##     69    185    200    204    218    243    244    252    256    260

Nous avons ici enlever les échantillons de moins de 1000 lectures. C’est
une manière d’éliminer les rares.

``` r
ps <- prune_samples(rowSums(otu_table(ps)) > 1000, ps)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
```

``` r
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")
evals <- out.pcoa.log$values[,1]
plot_ordination(pslog, out.pcoa.log, color = "age_binned",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->
Ce graphique a été réalisé par analyse PCoA (Principal COordinate
Analysis) en utilisant l’indice de dissimilarité de Bray-Curtis. Cela
permet notamment de visualiser la bêta-diversité en mesurant les
distances relatives entre communautés, à partir d’une matrice
d’observation. Pour revenir à l’analyse écologique, on peut ici
constater que l’âge et les portées ont un impact sur la communauté
microbienne. On peut en effet noter que plusieurs groupes se
distinguent.

``` r
library(phyloseq)
out.dpcoa.log <- ordinate(pslog, method = "DPCoA")
evals <- out.dpcoa.log$eig
plot_ordination(pslog, out.dpcoa.log, color = "age_binned", label= "SampleID",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->
Dans ce cas-là, il s’agit d’une analyse double de DPCoA. En effet, cette
méthode permet de représenter à la fois les échantillons et les
catégories taxonomiques. On relève que très peu de souris âgées sont
représentées ici. De plus, on peut ainsi supposer une liaison
score-genre taxonomique.

``` r
plot_ordination(pslog, out.dpcoa.log, type = "species", color = "Phylum") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->
Ce graphique incarne les taxons responsables des axes 1 et 2 du graphe
précedent. On peut ainsi dire que les échantillons responsables des
scores élevés de l’axe 2 contiennent plus de bactéroides. De même pour
l’axe 1 où les échantillons concernés contiennent un fort ensemble de
firmicutes.

``` r
out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCGAGCGTTATCCGGAATTATTGGGCGTAAAGAGGGAGCAGGCGGCAGCTAAGGTCTGCGGTGAAAGCCCGAAGCTAAACTTCGGTAAGCCGTGGAAACCGAGCAGCTAGAGTGCAGTAGAGGATCGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCAGTGGCGAAGGCGACGATCTGGGCTGCAACTGACGCTCAGTCCCGAAAGCGTGGGGAGC
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned",
                  shape = "family_relationship") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Binned Age", shape = "Litter")
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->
Ce graphique incarne les positions d’échantillons produites par un DPCoA
utilisant l’Unifrac pondéré. L’UniFrac pondéré tient compte de
l’abondance des organismes observés. Néanmoins, on relève entre autre
que l’interprétation de ce graphe est beaucoup moins précise que le
précédent. Toutefois, on remarque quand même une relation entre le
deuxième axe et l’âge des souris.

## Pourquoi les ‘plots’ d’ordination sont-ils si éloignés?

### Analyse PCA sur les rangs

``` r
abund <- otu_table(pslog) # nouvelle matrice représentant les abondances par leurs rangs
abund_ranks <- t(apply(abund, 1, rank)) # le microbe le plus petit dans un échantillon est mis en correspondance avec le rang 1
```

``` r
# tous les micro-organismes dont le rang est inférieur à un certain seuil sont fixés à 1
abund_ranks <- abund_ranks - 329
abund_ranks[abund_ranks < 1] <- 1 
```

Ces codes permettent donc de s’assurer qu’il n’y ai pas de différence de
rang en cas de taxons rares ou peu abondants.

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:Biostrings':
    ## 
    ##     collapse, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:XVector':
    ## 
    ##     slice

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(reshape2)
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

    ## Joining, by = c("Var1", "Var2")

``` r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

    ## Joining, by = c("Var1", "Var2")

``` r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

sample_ix <- sample(1:nrow(abund_df), 8)
ggplot(abund_df %>%
         filter(sample %in% abund_df$sample[sample_ix])) +
  geom_point(aes(x = abund, y = rank, col = sample),
             position = position_jitter(width = 0.2), size = 1.5) +
  labs(x = "Abundance", y = "Thresholded rank") +
  scale_color_brewer(palette = "Set2")
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-56-1.png)<!-- -->
On peut observer sur ce graphique le seuil de classement en fonction de
l’abondance. Il illustre les transformations effectuées dans les deux
blocs de codes ci-dessus. On peut deviner une relation presque
proportionnelle entre le rang et l’abondance. Toutefois, le début des
ensembles de points est marqué par un amas qui démontre bien que tous
les micro-organismes dont le rang est inférieur à un certain seuil sont
fixés à 1.

``` r
library(ade4)
```

    ## 
    ## Attaching package: 'ade4'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     score

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     score

``` r
ranks_pca <- dudi.pca(abund_ranks, scannf = F, nf = 3)
row_scores <- data.frame(li = ranks_pca$li,
                         SampleID = rownames(abund_ranks))
col_scores <- data.frame(co = ranks_pca$co,
                         seq = colnames(abund_ranks))
tax <- tax_table(ps) %>%
  data.frame(stringsAsFactors = FALSE)
tax$seq <- rownames(tax)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
tax$otu_id <- seq_len(ncol(otu_table(ps)))
row_scores <- row_scores %>%
  left_join(sample_data(pslog))
```

    ## Joining, by = "SampleID"

    ## Warning in class(x) <- c(setdiff(subclass, tibble_class), tibble_class): Setting
    ## class(x) to multiple strings ("tbl_df", "tbl", ...); result will no longer be an
    ## S4 object

``` r
col_scores <- col_scores %>%
  left_join(tax)
```

    ## Joining, by = "seq"

``` r
evals_prop <- 100 * (ranks_pca$eig / sum(ranks_pca$eig))
ggplot() +
  geom_point(data = row_scores, aes(x = li.Axis1, y = li.Axis2), shape = 2) +
  geom_point(data = col_scores, aes(x = 25 * co.Comp1, y = 25 * co.Comp2, col = Order),
             size = .3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  facet_grid(~ age_binned) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  coord_fixed(sqrt(ranks_pca$eig[2] / ranks_pca$eig[1])) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-58-1.png)<!-- -->
La figure ci-dessus incarne une analyse en composantes principales
(ACP). Ce type d’étude permet de caractériser puis visualiser un jeu de
données contenant des individus décrits par plusieurs variables. On peut
noter que les axes n’ont pas de dimension spéciale; ils représentent
plutôt des pourcentages de variance. On sait que plus les points sont
éloignés, plus les communautés seront différentes. De ce fait, on peut
dire que plus les souris vieillissent, plus les communautés microbiennes
se diversifient.

### Correspondance canonique

``` r
ps_ccpna <- ordinate(pslog, "CCA", formula = pslog ~ age_binned + family_relationship)
```

``` r
library(ggrepel)
ps_scores <- vegan::scores(ps_ccpna)
sites <- data.frame(ps_scores$sites)
sites$SampleID <- rownames(sites)
sites <- sites %>%
  left_join(sample_data(ps))
```

    ## Joining, by = "SampleID"

    ## Warning in class(x) <- c(setdiff(subclass, tibble_class), tibble_class): Setting
    ## class(x) to multiple strings ("tbl_df", "tbl", ...); result will no longer be an
    ## S4 object

``` r
species <- data.frame(ps_scores$species)
species$otu_id <- seq_along(colnames(otu_table(ps)))
species <- species %>%
  left_join(tax)
```

    ## Joining, by = "otu_id"

``` r
evals_prop <- 100 * ps_ccpna$CCA$eig[1:2] / sum(ps_ccpna$CA$eig)
ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
  geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 0.5) +
  geom_text_repel(data = species %>% filter(CCA2 < -2),
                    aes(x = CCA1, y = CCA2, label = otu_id),
            size = 1.5, segment.size = 0.1) +
  facet_grid(. ~ family_relationship) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
        y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  scale_color_brewer(palette = "Set2") +
  coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.45   ) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-60-1.png)<!-- -->
On peut relever que sur toutes les communautés précédemment remarquées,
seules quatres ont été annotées sur ces graphes.Cela permet de faciliter
notre lecture du bio-environnement. Comme précédemment, l’objectif est
ici de déterminer quelles communautés bactériennes sont les plus
importantes dans les différents types d’échantillons de souris. De cette
manière, les positions des échantillons sont déterminées par la
similarité des signatures des espèces et des caractéristiques
environnementales. Ces deux graphes nous permettent ici de noter la
présence d’une relation entre la portée des souris et leurs communautés
microbiennes. On relève aussi que les clostridies et ‘autres’ sont
vraiment peu abondants.

# Apprentissage supervisé

``` r
library(caret)
```

    ## Loading required package: lattice

``` r
library(ggplot2)
sample_data(pslog)$age2 <- cut(sample_data(pslog)$age, c(0, 100, 400))
dataMatrix <- data.frame(age = sample_data(pslog)$age2, otu_table(pslog))
# take 8 mice at random to be the training set, and the remaining 4 the test set
trainingMice <- sample(unique(sample_data(pslog)$host_subject_id), size = 8)
inTrain <- which(sample_data(pslog)$host_subject_id %in% trainingMice)
training <- dataMatrix[inTrain,]
testing <- dataMatrix[-inTrain,]
plsFit <- train(age ~ ., data = training,
                method = "pls", preProc = "center")
```

On a pu voir précédemment que les microbiomes évoluaient avec l’âge des
souris. On va ici s’atteler à prédire l’âge des souris à partir de la
composition de leurs microbiomes.

``` r
plsClasses <- predict(plsFit, newdata = testing)
table(plsClasses, testing$age)
```

    ##            
    ## plsClasses  (0,100] (100,400]
    ##   (0,100]        69         1
    ##   (100,400]       4        45

Ces résultats là traduisent un excellent travail de prédiction de l’âge.

``` r
library(randomForest) # forêts aléatoires comme autre exemple
```

    ## randomForest 4.6-14

    ## Type rfNews() to see new features/changes/bug fixes.

    ## 
    ## Attaching package: 'randomForest'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     margin

``` r
rfFit <- train(age ~ ., data = training, method = "rf",
               preProc = "center", proximity = TRUE)
rfClasses <- predict(rfFit, newdata = testing)
table(rfClasses, testing$age)
```

    ##            
    ## rfClasses   (0,100] (100,400]
    ##   (0,100]        71         1
    ##   (100,400]       2        45

``` r
library(vegan)
```

    ## Loading required package: permute

    ## This is vegan 2.5-6

    ## 
    ## Attaching package: 'vegan'

    ## The following object is masked from 'package:caret':
    ## 
    ##     tolerance

``` r
pls_biplot <- list("loadings" = loadings(plsFit$finalModel),
                   "scores" = scores(plsFit$finalModel))
class(pls_biplot$scores) <- "matrix"

pls_biplot$scores <- data.frame(sample_data(pslog)[inTrain, ],
                                pls_biplot$scores)

tax <- tax_table(ps)@.Data %>%
  data.frame(stringsAsFactors = FALSE)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
class(pls_biplot$loadings) <- "matrix"
pls_biplot$loadings <- data.frame(tax, pls_biplot$loadings)
```

``` r
ggplot() +
  geom_point(data = pls_biplot$scores,
             aes(x = Comp.1, y = Comp.2), shape = 2) +
  geom_point(data = pls_biplot$loadings,
             aes(x = 25 * Comp.1, y = 25 * Comp.2, col = Order),
             size = 0.3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Axis1", y = "Axis2", col = "Binned Age") +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  facet_grid( ~ age2) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-65-1.png)<!-- -->
On peut ici observer des graphiques conçus pour séparer les échantillons
par une variable. Il peut être interpréter de la même manière que les
précédents: c’est à dire que les communautés microbiennes varient avec
l’âge.

``` r
rf_prox <- cmdscale(1 - rfFit$finalModel$proximity) %>%
  data.frame(sample_data(pslog)[inTrain, ])

ggplot(rf_prox) +
  geom_point(aes(x = X1, y = X2, col = age_binned),
             size = 1, alpha = 0.7) +
  scale_color_manual(values = c("#A66EB8", "#238DB5", "#748B4F")) +
  guides(col = guide_legend(override.aes = list(size = 4))) +
  labs(col = "Binned Age", x = "Axis1", y = "Axis2")
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-66-1.png)<!-- -->
Le modèle de forêt aléatoire détermine ici une distance entre les
échantillons, qui peut être entrée dans le PCoA pour produire un ‘plot’
de proximité. Autrement dit, ce graphique représente une courbe de
proximité aléatoire de la forêt. On sait que dans le cas où les
échantillons sont fréquemment rencontrés, il y aura une faible distance
entre eux sur le graphique. On peut donc dire que les classes d’arbres
sont clairement séparées par leur âge.

``` r
as.vector(tax_table(ps)[which.max(importance(rfFit$finalModel)), c("Family", "Genus")])
```

    ## [1] "Lachnospiraceae" "Roseburia"

Nous avons d’abord identifer le micro-organisme ayant le plus
d’influence dans la forêt aléatoire. Il s’agit de la famille des
Lachnospiraceae, du genre Roseburia.

``` r
impOtu <- as.vector(otu_table(pslog)[,which.max(importance(rfFit$finalModel))])
maxImpDF <- data.frame(sample_data(pslog), abund = impOtu)
ggplot(maxImpDF) +   geom_histogram(aes(x = abund)) +
  facet_grid(age2 ~ .) +
  labs(x = "Abundance of discriminative bacteria", y = "Number of samples")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-68-1.png)<!-- -->
Cet histogramme permet de visualiser l’abondance de ce micro-organisme.
On peut voir qu’il est uniformément peu présent de 0 à 100 jours (à
l’exception des tous premiers). Au contraire, son abondance est plus
élevée de 100 à 400 jours.

# Analyses basées sur des graphiques

## Créer et tracer des graphiques

Ci-dessous, nous avons créer un réseau grâce à l’indice de similarité de
Jaccard. Ce dernier ne tenant pas compte de l’abondance, une seule
présence contera autant que plusieurs.

``` r
library("phyloseqGraphTest")
library("igraph")
```

    ## 
    ## Attaching package: 'igraph'

    ## The following object is masked from 'package:vegan':
    ## 
    ##     diversity

    ## The following object is masked from 'package:permute':
    ## 
    ##     permute

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     as_data_frame, groups, union

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     union

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     union

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     normalize, path, union

    ## The following objects are masked from 'package:ape':
    ## 
    ##     edges, mst, ring

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
library("ggnetwork")
net <- make_network(ps, max.dist=0.35)
sampledata <- data.frame(sample_data(ps))
V(net)$id <- sampledata[names(V(net)), "host_subject_id"]
V(net)$litter <- sampledata[names(V(net)), "family_relationship"]

net_graph <- ggnetwork(net)

ggplot(net_graph, aes(x = x, y = y, xend = xend, yend = yend), layout = "fruchtermanreingold") +
  geom_edges(color = "darkgray") +
  geom_nodes(aes(color = id, shape = litter),  size = 3 ) +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        legend.key.height = unit(0.5,"line")) +
  guides(col = guide_legend(override.aes = list(size = .5)))
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-69-1.png)<!-- -->
Ce réseau a donc été créé par le seuillage de la matrice de
dissimilarité (ou similarité) Jaccard. Les couleurs incarnent la souris
d’où provient l’échantillon et la forme représente la portée dans
laquelle se trouvait cette même souris. On sait que plus il y a de
points connectés, plus les communautés partagent des micro-organismes
identiques. Nous pouvons voir qu’il y a un regroupement des échantillons
à la fois par souris et par portée. Néanmoins, mis à part le réseau
principal, on note la présence de plusieurs duos à quatuors,
représentant les variants. On peut également remarquer qu’il est
possible que des souris de portées différents soient liées par leur
microbiote, même si ce n’est pas majoritairement le cas.

## Tests à deux échantillons basés sur des graphiques

### Arbre à Portée Minimale

Un arbre à portée minimale est basé premièrement sur les distances entre
les échantillons, puis sur le comptage du nombre d’arêtes de l’arbre qui
se trouvaient entre les échantillons de différents groupes.

``` r
gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
                      distance = "jaccard", type = "mst")
gt$pval
```

    ## [1] 0.004

``` r
plotNet1=plot_test_network(gt) + theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
plotPerm1=plot_permutations(gt)
grid.arrange(ncol = 2,  plotNet1, plotPerm1)
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-71-1.png)<!-- -->
Ce histogramme a été obtenu à partir de l’arbre à portée minimale défini
grâce à l’indice de Jaccard. On peut remarquer grâce à cet arbre que les
échantillons se regroupent bien par portée. Un regroupement mixte est
marqué par les deux couleurs différentes et par la ligne en pointillés
les reliant. On peut noter sur l’histogramme un pic au niveau de 245
arrêtes pures. Il signifie qu’avec une soixantaine de comptages, le pic
de pureté des communautés est atteint.

### Les plus proches voisins

``` r
gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
                      distance = "jaccard", type = "knn", knn = 1)
```

``` r
plotNet2=plot_test_network(gt) + theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
plotPerm2=plot_permutations(gt)
grid.arrange(ncol = 2,  plotNet2, plotPerm2)
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-73-1.png)<!-- -->
Sur le graphique, il est fort probable que lorsque deux échantillons
soient reliés, ils soient de la même portée. Il s’agit d’un réseau de
voisins proches (avec k=1). De la même manière que le graphique
précédent, un regroupement mixte est marqué par les deux couleurs
différentes et par la ligne en pointillés les reliant. Néanmoins, avec
k=1 dans ce cas, on peut considéré la diminution des possibilités de
combinaison ou de rapprochement par exemple. Par ailleurs, l’histogramme
montre un pic de comptage à environ 50, ce qui correspond à à peu près
205 arrêtes pures. Autrement dit, il faut 50 comptages environ pour
avoir un type d’arrêtes pures à environ 205 sur le graphique de gauche.

## Modélisation linéaire

``` r
library("nlme")
```

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     collapse

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     collapse

``` r
library("reshape2")
ps_alpha_div <- estimate_richness(ps, split = TRUE, measure = "Shannon")
ps_alpha_div$SampleID <- rownames(ps_alpha_div) %>%
  as.factor()
ps_samp <- sample_data(ps) %>%
  unclass() %>%
  data.frame() %>%
  left_join(ps_alpha_div, by = "SampleID") %>%
  melt(measure.vars = "Shannon",
       variable.name = "diversity_measure",
       value.name = "alpha_diversity")

# réorganiser la facette de la diversité de la plus basse à la plus haute
diversity_means <- ps_samp %>%
  group_by(host_subject_id) %>%
  summarise(mean_div = mean(alpha_diversity)) %>%
  arrange(mean_div)
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
ps_samp$host_subject_id <- factor(ps_samp$host_subject_id)
#                                  diversity_means$host_subject_id)
```

``` r
alpha_div_model <- lme(fixed = alpha_diversity ~ age_binned, data = ps_samp,
                       random = ~ 1 | host_subject_id)
```

``` r
new_data <- expand.grid(host_subject_id = levels(ps_samp$host_subject_id),
                        age_binned = levels(ps_samp$age_binned))
new_data$pred <- predict(alpha_div_model, newdata = new_data)
X <- model.matrix(eval(eval(alpha_div_model$call$fixed)[-2]),
                  new_data[-ncol(new_data)])
pred_var_fixed <- diag(X %*% alpha_div_model$varFix %*% t(X))
new_data$pred_var <- pred_var_fixed + alpha_div_model$sigma ^ 2
```

``` r
# fitted values, with error bars
ggplot(ps_samp %>% left_join(new_data)) +
  geom_errorbar(aes(x = age_binned, ymin = pred - 2 * sqrt(pred_var),
                    ymax = pred + 2 * sqrt(pred_var)),
                col = "#858585", size = .1) +
  geom_point(aes(x = age_binned, y = alpha_diversity,
                 col = family_relationship), size = 0.8) +
  facet_wrap(~host_subject_id) +
  scale_y_continuous(limits = c(2.4, 4.6), breaks = seq(0, 5, .5)) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Binned Age", y = "Shannon Diversity", color = "Litter") +
  guides(col = guide_legend(override.aes = list(size = 4))) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)),
        axis.text.x = element_text(angle = -90, size = 6),
        axis.text.y = element_text(size = 6))
```

    ## Joining, by = c("host_subject_id", "age_binned")

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-77-1.png)<!-- -->
Cet ensemble de graphes représente l’alpha-diversité (par l’indice de
Shannon) en fonction de l’âge. L’utilisation de cet indice permet de
refléter au mieux l’abondance des communautés. Ces graphiques permettent
d’étudier la relation entre la diversité des communautés microbiennes de
souris, les variables d’âge et celles de portée. On peut dans un premier
temps remarquer que les jeunes souris présentent une alpha-diversité de
Shannon sensiblement plus faible. Ensuite, on peut confirmer le fait que
les communautés microbiennes sont étroitement liées à la varialbe de
portée. Finalement, on relève que chez les souris agées il y a
extrêmement peu, voire pas du tout, de diversité de Shannon, soit de
diversité de communauté microbienne.

## Tests multiples hiérarchisés

``` r
library("reshape2")
library("DESeq2")
```

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     sampleNames

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following object is masked from 'package:igraph':
    ## 
    ##     simplify

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply, rowsum

``` r
#New version of DESeq2 needs special levels
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship = gsub(" ", "", sample_data(ps)$family_relationship)
ps_dds <- phyloseq_to_deseq2(ps, design = ~ age_binned + family_relationship)
```

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` r
# geometric mean, set to zero when all coordinates are zero
geo_mean_protected <- function(x) {
  if (all(x == 0)) {
    return (0)
  }
  exp(mean(log(x[x != 0])))
}

geoMeans <- apply(counts(ps_dds), 1, geo_mean_protected)
ps_dds <- estimateSizeFactors(ps_dds, geoMeans = geoMeans)
ps_dds <- estimateDispersions(ps_dds)
```

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

``` r
abund <- getVarianceStabilizedData(ps_dds)
```

``` r
short_names <- substr(rownames(abund), 1, 5)%>%
  make.names(unique = TRUE)
rownames(abund) <- short_names
```

``` r
abund_sums <- rbind(data.frame(sum = colSums(abund),
                               sample = colnames(abund),
                               type = "DESeq2"),
                    data.frame(sum = rowSums(otu_table(pslog)),
                               sample = rownames(otu_table(pslog)),
                               type = "log(1 + x)"))

ggplot(abund_sums) +
  geom_histogram(aes(x = sum), binwidth = 20) +
  facet_grid(type ~ .) +
  xlab("Total abundance within sample")
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-80-1.png)<!-- -->
Ces deux histogrammes montrent deux méthodes différentes de modification
de population : par DeSeq2 et logarithme. On peut tout d’abord relever
que les ensembles microbiens ont la même allure. Toutefois, on peut
aussi noter que le graphique par rapport à DeSeq2 est plus étalé sur les
faibles valeurs d’abondance totale, contrairement à celui en logarithme
qui s’étale plutôt sur les valeurs élevées. Il semble donc que la
méthode DeSeq2 soit plus efficace. En effet, avec un comptage d’à peu
près 35, l’abondance totale est encore aux environs des 600, ce qui ne
sera pas le cas ici avec le logarithme.

``` r
library("structSSI")
el <- phy_tree(pslog)$edge
el0 <- el
el0 <- el0[nrow(el):1, ]
el_names <- c(short_names, seq_len(phy_tree(pslog)$Nnode))
el[, 1] <- el_names[el0[, 1]]
el[, 2] <- el_names[as.numeric(el0[, 2])]
unadj_p <- treePValues(el, abund, sample_data(pslog)$age_binned)
```

``` r
hfdr_res <- hFDR.adjust(unadj_p, el, .75)
summary(hfdr_res)
```

    ## Number of hypotheses: 764 
    ## Number of tree discoveries: 579 
    ## Estimated tree FDR: 1 
    ## Number of tip discoveries: 280 
    ## Estimated tips FDR: 1 
    ## 
    ##  hFDR adjusted p-values: 
    ##                 unadjp         adjp adj.significance
    ## GCAAG.95  1.861873e-82 3.723745e-82              ***
    ## GCAAG.70  1.131975e-75 2.263950e-75              ***
    ## GCAAG.187 5.148758e-59 1.029752e-58              ***
    ## GCAAG.251 3.519276e-50 7.038553e-50              ***
    ## GCAAG.148 1.274481e-49 2.548962e-49              ***
    ## GCAAG.30  9.925218e-49 1.985044e-48              ***
    ## GCGAG.76  1.722591e-46 3.445183e-46              ***
    ## GCAAG.167 6.249050e-43 1.249810e-42              ***
    ## 255       8.785479e-40 1.757096e-39              ***
    ## GCAAG.64  2.727610e-36 5.455219e-36              ***
    ## [only 10 most significant hypotheses shown] 
    ## --- 
    ## Signif. codes:  0 '***' 0.015 '**' 0.15 '*' 0.75 '.' 1.5 '-' 1

``` r
# partie à ne pas éxécuter
plot(hfdr_res, height = 5000) # s'ouvre dans un navigateur
```

Normalement, les lignes de codes et résultats ci-dessus permettent
d’obtenir un graphique avec de nombreuses bactéries (GCAAG.95,
GCAAG.70, etc.), d’abondance différente. Il ne peut toutefois ici pas
s’ouvrir dans un navigateur.

``` r
tax <- tax_table(pslog)[, c("Family", "Genus")] %>%
  data.frame()
tax$seq <- short_names
```

``` r
options(digits=3)
hfdr_res@p.vals$seq <- rownames(hfdr_res@p.vals)
tax %>%
  left_join(hfdr_res@p.vals) %>%
  arrange(adjp) %>% head(10)
```

    ## Joining, by = "seq"

    ##             Family            Genus       seq   unadjp     adjp
    ## 1  Lachnospiraceae             <NA>  GCAAG.95 1.86e-82 3.72e-82
    ## 2  Lachnospiraceae        Roseburia  GCAAG.70 1.13e-75 2.26e-75
    ## 3  Lachnospiraceae Clostridium_XlVa GCAAG.187 5.15e-59 1.03e-58
    ## 4  Lachnospiraceae             <NA> GCAAG.251 3.52e-50 7.04e-50
    ## 5  Lachnospiraceae Clostridium_XlVa GCAAG.148 1.27e-49 2.55e-49
    ## 6  Lachnospiraceae             <NA>  GCAAG.30 9.93e-49 1.99e-48
    ## 7  Ruminococcaceae     Ruminococcus  GCGAG.76 1.72e-46 3.45e-46
    ## 8  Lachnospiraceae Clostridium_XlVa GCAAG.167 6.25e-43 1.25e-42
    ## 9  Lachnospiraceae        Roseburia  GCAAG.64 2.73e-36 5.46e-36
    ## 10            <NA>             <NA>   GCAAG.1 5.22e-35 1.04e-34
    ##    adj.significance
    ## 1               ***
    ## 2               ***
    ## 3               ***
    ## 4               ***
    ## 5               ***
    ## 6               ***
    ## 7               ***
    ## 8               ***
    ## 9               ***
    ## 10              ***

On peut voir d’après les données ci-dessus que les bactéries le plus
souvent associées sont des Lachnospiraceae. Ces résultats correspondent
donc à ceux de la forêt aléatoire. En effet, on note la présence entre
autre de Roseburia et Clostridium.

# ‘Multitable techniques’

Tout d’abord, les codes ci-dessous récupèrent puis filtrent un nouvel
ensemble de données.

``` r
metab <- read.csv("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/metabolites.csv",row.names = 1)
microbe_connect <-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/microbe.rda")
load(microbe_connect)
microbe
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 20609 taxa and 12 samples ]
    ## tax_table()   Taxonomy Table:    [ 20609 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 20609 tips and 20607 internal nodes ]

Il s’agit là de filtrer les données nulles de nombreux échantillons.

``` r
library("genefilter")
```

    ## 
    ## Attaching package: 'genefilter'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     rowSds, rowVars

``` r
keep_ix <- rowSums(metab == 0) <= 3
metab <- metab[keep_ix, ]
microbe <- prune_taxa(taxa_sums(microbe) > 4, microbe)
microbe <- filter_taxa(microbe, filterfun(kOverA(3, 2)), TRUE)
metab <- log(1 + metab, base = 10)
X <- otu_table(microbe)
X[X > 50] <- 50
dim(X)
```

    ## [1] 174  12

``` r
dim(metab)
```

    ## [1] 405  12

Ensuite, une Analyse Canonique par Corrélation (CCA) permet de comparer
les données des deux tableaux crées X et metab (contenant chacun 12
colonnes). Ce type d’analyse compare deux groupes de variables pour
savoir s’ils décrivent un même phénomène, auquel cas on pourra se passer
de l’un des deux (procédure de sélection).

``` r
library(PMA)
cca_res <- CCA(t(X),  t(metab), penaltyx = .15, penaltyz = .15)
```

    ## 123456789101112131415

``` r
cca_res
```

    ## Call: CCA(x = t(X), z = t(metab), penaltyx = 0.15, penaltyz = 0.15)
    ## 
    ## 
    ## Num non-zeros u's:  5 
    ## Num non-zeros v's:  15 
    ## Type of x:  standard 
    ## Type of z:  standard 
    ## Penalty for x: L1 bound is  0.15 
    ## Penalty for z: L1 bound is  0.15 
    ## Cor(Xu,Zv):  0.974

Les deux tableaux analysés sont corrélés à 97,4%. Ils sont donc
approximativement similaires.

``` r
combined <- cbind(t(X[cca_res$u != 0, ]),
                  t(metab[cca_res$v != 0, ]))
pca_res <- dudi.pca(combined, scannf = F, nf = 3)
```

``` r
genotype <- substr(rownames(pca_res$li), 1, 2)
sample_type <- substr(rownames(pca_res$l1), 3, 4)
feature_type <- grepl("\\.", colnames(combined))
feature_type <- ifelse(feature_type, "Metabolite", "OTU")
sample_info <- data.frame(pca_res$li, genotype, sample_type)
feature_info <- data.frame(pca_res$c1,
                           feature = substr(colnames(combined), 1, 6))
```

``` r
ggplot() +  geom_point(data = sample_info,
            aes(x = Axis1, y = Axis2, col = sample_type, shape = genotype), size = 3) + 
  geom_label_repel(data = feature_info,
                   aes(x = 5.5 * CS1, y = 5.5 * CS2, label = feature, fill = feature_type),
                   size = 2, segment.size = 0.3,
                   label.padding = unit(0.1, "lines"), label.size = 0) +
  geom_point(data = feature_info,
             aes(x = 5.5 * CS1, y = 5.5 * CS2, fill = feature_type),
             size = 1, shape = 23, col = "#383838") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#a6d854", "#e78ac3")) +
  guides(fill = guide_legend(override.aes = list(shape = 32, size = 0))) +
  coord_fixed(sqrt(pca_res$eig[2] / pca_res$eig[2])) +
  labs(x = sprintf("Axis1 [%s%% Variance]",
                   100 * round(pca_res$eig[1] / sum(pca_res$eig), 2)),
       y = sprintf("Axis2 [%s%% Variance]",
                   100 * round(pca_res$eig[2] / sum(pca_res$eig), 2)),
       fill = "Feature Type", col = "Sample Type")
```

![](03_data_analysis_CC1_files/figure-gfm/unnamed-chunk-93-1.png)<!-- -->
Ce graphique par PCA montre trois caractéristiques : le type
d’échantillon, métabolite ou OTU, et enfin le génotype (knockout ou
type sauvage). On peut voir parmi ces trois caractéristiques que seul le
type d’échantillon (différents régimes alimentaires) permet de
différencier distinctement les communautés bactériennes. On peut donc
dire que le régime alimentaire influe fortement sur les communautés
microbiennes de l’organisme.
