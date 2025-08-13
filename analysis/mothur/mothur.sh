#################################################################
###                   Krys Kibler 2024-03-12                  ###
### Purpose: To process 16S DNA from LView test-run samples   ###
### Resource: https://mothur.org/wiki/miseq_sop/              ###
#################################################################

# Done in interactive condor on glbrc server

# Pipelines
MOTHUR=/opt/bifxapps/mothur-1.48.0

# Code
cd /home/glbrc.org/kjkibler/LAKEVIEW-microbe/data/sequences/trial-run_LakeView16S

# activate mothur interactive
# copy lines after `$MOTHUR/mothur >` into terminal for running code
$MOTHUR/mothur

# Which fastq files go together
# output: stability.files
$MOTHUR/mothur > make.file(inputdir=., type=fastq, prefix=stability)


#################################################################
# Reducing sequencing and PCR errors

# Make 16S "contigs"
$MOTHUR/mothur > make.contigs(file=stability.files)

#Group count:
#Lakeview1	41897
#Lakeview2	37722
#Lakeview3	40289
#Lakeview4	66489
#Lakeview5	47088
#Lakeview6	41634
#Lakeview7	44451
#Lakeview8	42397

#Total of all groups is 361967

#It took 104 secs to process 361967 sequences.

#Output File Names:
#stability.trim.contigs.fasta
#stability.scrap.contigs.fasta
#stability.contigs_report
#stability.contigs.count_table

# check out sequences
$MOTHUR/mothur > summary.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table)

#		Start	End	NBases	Ambigs	Polymer	NumSeqs
#Minimum:	1	248	248	0	3	1
#2.5%-tile:	1	253	253	0	4	9050
#25%-tile:	1	253	253	0	4	90492
#Median: 	1	253	253	0	4	180984
#75%-tile:	1	253	253	0	4	271476
#97.5%-tile:	1	254	254	1	6	352918
#Maximum:	1	500	500	70	228	361967
#Mean:	1	253	253	0	4
# of unique seqs:	361967
#total # of seqs:	361967

#It took 9 secs to summarize 361967 sequences.

#Output File Names:
#stability.trim.contigs.summary


# remove sus "contigs"
$MOTHUR/mothur > screen.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table, maxambig=0, maxlength=275, maxhomop=8)

# check out sequences
$MOTHUR/mothur > summary.seqs(count=current)

#		Start	End	NBases	Ambigs	Polymer	NumSeqs
#Minimum:	1	250	250	0	3	1
#2.5%-tile:	1	253	253	0	4	8560
#25%-tile:	1	253	253	0	4	85593
#Median: 	1	253	253	0	4	171186
#75%-tile:	1	253	253	0	4	256779
#97.5%-tile:	1	254	254	0	6	333812
#Maximum:	1	275	275	0	8	342371
#Mean:	1	253	253	0	4
# of unique seqs:	342371
#total # of seqs:	342371

#It took 9 secs to summarize 342371 sequences.

#Output File Names:
#stability.trim.contigs.good.summary

#################################################################
# Processing improved sequences

# curate only unique "contigs"
$MOTHUR/mothur > unique.seqs(fasta=stability.trim.contigs.good.fasta, count=stability.contigs.good.count_table)

# summary stats
$MOTHUR/mothur > summary.seqs(count=stability.trim.contigs.good.count_table)

# of unique seqs:	47403
#total # of seqs:	342371


# customize alignment reference based on region amplified (use coordingates; start=11895, end=25318 == V4 region)
# https://mothur.org/blog/2016/Customization-for-your-region/ -- blog post to help figure out bp coordinates
$MOTHUR/mothur > pcr.seqs(fasta=silva.nr_v132.align, start=11895, end=25318, keepdots=F)

# rename reference alignment
$MOTHUR/mothur > rename.file(input=silva.nr_v132.pcr.align, new=silva.v4.fasta)

# align sequences based on made reference
$MOTHUR/mothur > align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=silva.v4.fasta)

# generate summary table to check alignments
$MOTHUR/mothur > summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table)

#		Start	End	NBases	Ambigs	Polymer	NumSeqs
#Minimum:	1	829	2	0	1	1
#2.5%-tile:	1968	11550	253	0	4	8560
#25%-tile:	1968	11550	253	0	4	85593
#Median: 	1968	11550	253	0	4	171186
#75%-tile:	1968	11550	253	0	4	256779
#97.5%-tile:	1968	11550	254	0	6	333812
#Maximum:	13422	13424	275	0	8	342371
#Mean:	1968	11549	253	0	4
# of unique seqs:	47403
#total # of seqs:	342371

#It took 6 secs to summarize 342371 sequences.


# re-sreen aligned "contigs" to only get the ones that are good
$MOTHUR/mothur > screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, start=1968, end=11550)
#It took 6 secs to screen 47403 sequences, removed 301.

$MOTHUR/mothur > summary.seqs(fasta=current, count=current)

#		Start	End	NBases	Ambigs	Polymer	NumSeqs
#Minimum:	1	11550	250	0	3	1
#2.5%-tile:	1968	11550	253	0	4	8548
#25%-tile:	1968	11550	253	0	4	85479
#Median: 	1968	11550	253	0	4	170957
#75%-tile:	1968	11550	253	0	4	256435
#97.5%-tile:	1968	11550	254	0	6	333366
#Maximum:	1968	13424	272	0	8	341913
#Mean:	1967	11550	253	0	4
# of unique seqs:	47102
#total # of seqs:	341913

#It took 6 secs to summarize 341913 sequences.

# trims the ends
$MOTHUR/mothur > filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)

#Length of filtered alignment: 495
#Number of columns removed: 12929
#Length of the original alignment: 13424
#Number of sequences used to construct filter: 47102

# re-curate for unique reads
$MOTHUR/mothur > unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)

# clustering
$MOTHUR/mothur > pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=2)

# check out
$MOTHUR/mothur > summary.seqs(fasta=current, count=current)

#		Start	End	NBases	Ambigs	Polymer	NumSeqs
#Minimum:	1	494	223	0	3	1
#2.5%-tile:	1	495	253	0	4	8548
#25%-tile:	1	495	253	0	4	85479
#Median: 	1	495	253	0	4	170957
#75%-tile:	1	495	253	0	4	256435
#97.5%-tile:	1	495	254	0	6	333366
#Maximum:	2	495	268	0	8	341913
#Mean:	1	494	253	0	4
# of unique seqs:	18217
#total # of seqs:	341913


# removing sequences with suspected chimearas
$MOTHUR/mothur > chimera.vsearch(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)

# check to see how many got removed
$MOTHUR/mothur > summary.seqs(fasta=current, count=current)

#		Start	End	NBases	Ambigs	Polymer	NumSeqs
#Minimum:	1	494	223	0	3	1
#2.5%-tile:	1	495	253	0	4	8180
#25%-tile:	1	495	253	0	4	81800
#Median: 	1	495	253	0	4	163599
#75%-tile:	1	495	253	0	4	245398
#97.5%-tile:	1	495	254	0	6	319017
#Maximum:	2	495	268	0	8	327196
#Mean:	1	494	253	0	4
# of unique seqs:	8907
#total # of seqs:	327196


# CLASSIFICATION
$MOTHUR/mothur > classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.count_table, reference=trainset18_062020.pds.fasta, taxonomy=trainset18_062020.pds.tax)
#or
$MOTHUR/mothur > classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.count_table, reference=silva.nr_v132.align, taxonomy=silva.nr_v132.tax)



# remove undesirable origin "contigs" (ie. eukaryotic and archae stuff)
$MOTHUR/mothur > remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
$MOTHUR/mothur > remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.nr_v132.wang.taxonomy, taxon=unknown)

# FINAL OUTPUT OF CLASSIFIED SEQUENCES
$MOTHUR/mothur > summary.tax(taxonomy=current, count=current)
# output == stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.nr_v132.wang.pick.tax.summary


#################################################################
# Generate OTUs

$MOTHUR/mothur > rename.file(fasta=current, count=current, taxonomy=current, prefix=final)

# (1) Good for small datasets; the "traditional" way
$MOTHUR/mothur > dist.seqs(fasta=final.fasta, cutoff=0.03)
$MOTHUR/mothur > cluster(column=final.dist, count=final.count_table)

# or (2) "Non-traditional" way
$MOTHUR/mothur > cluster.split(fasta=final.fasta, count=final.count_table, taxonomy=final.taxonomy, taxlevel=4, cutoff=0.03)

# How many "contigs" are in each OTU
$MOTHUR/mothur > make.shared(list=final.opti_mcc.list, count=final.count_table, label=0.03)

# Consensus classification
# output: final.opti_mcc.0.03.cons.taxonomy, final.opti_mcc.0.03.cons.tax.summary
$MOTHUR/mothur > classify.otu(list=final.opti_mcc.list, count=final.count_table, taxonomy=final.taxonomy, label=0.03)


#################################################################
# Phylogenetics
$MOTHUR/mothur > dist.seqs(fasta=final.fasta, output=lt)
$MOTHUR/mothur > clearcut(phylip=final.phylip.dist)

# output: final.phylip.dist

#################################################################
# Analysis

# count groups to "rarify" samples
$MOTHUR/mothur > count.groups(shared=final.opti_mcc.shared)

#Lakeview1 contains 37614.
#Lakeview2 contains 34655.
#Lakeview3 contains 37521.
#Lakeview4 contains 56026.
#Lakeview5 contains 43530.
#Lakeview6 contains 38255.
#Lakeview7 contains 40629.
#Lakeview8 contains 37851.

#Size of smallest group: 34655.

#Total seqs: 326081.

#Output File Names:
#final.opti_mcc.count.summary

mothur > sub.sample(shared=final.opti_mcc.shared, size=34655)

# OTU based Analysis

# Alpha diversity
mothur > rarefaction.single(shared=final.opti_mcc.shared, calc=sobs, freq=100)

#Output File Names:
#final.opti_mcc.groups.rarefaction

#number of sequences, the sample coverage, the number of observed OTUs, and the Inverse Simpson diversity estimate
mothur > summary.single(shared=final.opti_mcc.shared, calc=nseqs-coverage-sobs-invsimpson, subsample=T)

#Output File Names:
#final.opti_mcc.groups.ave-std.summary

# Beta diversity
mothur > dist.shared(shared=final.opti_mcc.shared, calc=braycurtis-jclass, subsample=t)

#Output File Names:
#final.opti_mcc.braycurtis.0.03.lt.ave.dist
#final.opti_mcc.braycurtis.0.03.lt.std.dist
#final.opti_mcc.jclass.0.03.lt.ave.dist
#final.opti_mcc.jclass.0.03.lt.std.dist

# ordination analysis
mothur > pcoa(phylip=final.opti_mcc.braycurtis.0.03.lt.ave.dist)

#Output File Names:
#final.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes
#final.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings

mothur > nmds(phylip=final.opti_mcc.braycurtis.0.03.lt.ave.dist)

#Number of dimensions:	2
#Lowest stress :	0.083167
#R-squared for configuration:	0.984752

#Output File Names:
#final.opti_mcc.braycurtis.0.03.lt.ave.nmds.iters
#final.opti_mcc.braycurtis.0.03.lt.ave.nmds.stress
#final.opti_mcc.braycurtis.0.03.lt.ave.nmds.axes

mothur > nmds(phylip=final.opti_mcc.braycurtis.0.03.lt.ave.dist, mindim=3, maxdim=3)


# anova
# need separate file for "treatments"
mothur > amova(phylip=final.opti_mcc.braycurtis.0.03.lt.ave.dist, design=lakeview.design.txt)
mothur > homova(phylip=final.opti_mcc.braycurtis.0.03.lt.ave.dist, design=lakeview.design.txt)

# which OTUs shift the axises
mothur > corr.axes(axes=final.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes, shared=final.opti_mcc.0.03.subsample.shared, method=spearman, numaxes=3)

#Output File Names:
#final.opti_mcc.0.03.subsample.spearman.corr.axes


# which metadata shifts the axises
# need metadata
mothur > corr.axes(axes=final.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes, metadata=mouse.dpw.metadata, method=spearman, numaxes=3)

# can we separate out different communities
mothur > get.communitytype(shared=final.opti_mcc.0.03.subsample.shared)



# Population level Analysis
mothur > metastats(shared=final.opti_mcc.0.03.subsample.shared, design=mouse.time.design)



# lefse
mothur > lefse(shared=final.opti_mcc.0.03.subsample.shared, design=lakeview.design.txt)

#No features with significant differences between the classes.

#Output File Names:
#final.opti_mcc.0.03.subsample.0.03.lefse_summary
