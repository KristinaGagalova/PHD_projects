# Use only the transcripts that are in common between the 2 target sets (PG29 and WS-Q903)

Extract the overlapping transcripts
```
while read line; do grep $line ParseSam/Out.results.WS77111.extra.top_sam.out; done < OverlapTranscripts1.txt > OverlapSucceededWS77111extraAll_sam.txt
while read line; do grep $line ParseSam/Out.results.Q903.top_sam.out; done < OverlapTranscripts1.txt > OverlapSucceededQ903_sam.txt
while read line; do grep $line ../../PreliminaryResults12May/data/ParseSam/Out.results.PG29.top_sam.out; done < OverlapTranscripts1.txt > OverlapTranscriptsPG29_sam.txt
```
