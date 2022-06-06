## Trim and process the reads to avoid sequencing errors

The directory describes a very light trimming.

1) Process the reads with [PRINseq](http://prinseq.sourceforge.net/manual.html). Filter reads with low complexity, reads containing a high "N" content. Trim the reads mildly for Phred score lower than 20.

2) Because beind paired end reads, some reads will be singletons after the trimming/filtering. Process those and sort them

3) Add fak epairs to the sinletons, this infor may be useful for BBT maker.

4) CheckIntegrity - calculates the number of lines, the number of reads and the lines with only "+" per fastq file

More details in [Jira](https://www.bcgsc.ca/jira/browse/BTL-908?focusedCommentId=823967&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-823967)
