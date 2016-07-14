# clj-blast

A parser for BLAST XML files.

## Usage

Import from Clojars:

```clojure
[clj-blast "0.2.7"]
```

Use in your namespace:

```clojure
(require '[clj-blast.core :as bl])
```

Open a reader on a blast output file and call 'iteration-seq' on the
reader to get a lazy list of iterations. To get hits call 'hit-seq' on
iterations to get maps representing the hits. HSPs can be accessed in
the :hsps field and alignments can be accessed in the :alignment field
of individual HSPs. :evalue or :bit-score keywords can be used to
filter hits returned by 'hit-seq' for signficance by ensuring at least
one HSP has a score higher (or lower for evalue) than the specified
scoring parameter.

Currently does not support Blasts XML2 format.

```clojure
user> (with-open [r (reader tf)]
                  (->> (iteration-seq r)
                       (mapcat hit-seq)
                       first))
{:Hit_id "sp|Q9BPA9|CO26_CONTE", :Hit_len 70, :Hit_accession "Q9BPA9",
 :Hit_def "Conotoxin 6 OS=Conus textile PE=1 SV=1", :Hit_num 1, :hsps
 ({:Hsp_query-from 9, :Hsp_score 50.0, :Hsp_midline "+P S+ APCC+  T + R N",
 :Hsp_density nil, :alignment "9   SPSSTRAPCCNSKTPATRVN  28\n    +P S+ APCC+ 
 T + R N  \n48  APCSSGAPCCDWWTCSARTN  67\n", :Hsp_identity 9, :Hsp_hit-from 48,
 :Hsp_hit-to 67, :Hsp_hit-frame 0, :Hsp_pattern-to nil, :Hsp_positive 13,
 :Hsp_align-len 20, :Hsp_query-frame 0, :Hsp_qseq "SPSSTRAPCCNSKTPATRVN",
 :Hsp_num 1, :Hsp_pattern-from nil, :Hsp_bit-score 23.8682, :Hsp_query-to 28,
 :Hsp_gaps 0, :Hsp_evalue 1.36307, :Hsp_hseq "APCSSGAPCCDWWTCSARTN"}),
 :query-accession "Query_1"}
user> (with-open [r (reader tf)]
                  (->> (iteration-seq r)
                       (mapcat #(hit-seq % :bit-score 40))
                       first))
{:Hit_id "sp|J3SDX8|LICH_CROAD", :Hit_len 400, :Hit_accession "J3SDX8", ...
user>
```

Blasts can be performed using `blast` and `blast-file`.

Sequences can be retrieved using `retrieve-sequence` and
`blastdb->file`:

```clojure
user> (with-open [r (reader sp)]
        (blastdb->file (take 10000 (map :accession (fasta-seq r)))
	               "/path/blastdb"
		       "/path/outfile"
"/path/outfile"
user>
user> (retrieve-sequence ["comp0_c0_seq1"] "/path/blastdb" "nucl")
({:accession "lcl|comp0_c0_seq1", :description "len=203 path=[521:0-202]",
 :sequence "GCGCATT..."})
user>
```

Note that `retrieve-sequence` is **not** lazy.

Blast databases can be created using `create-blastdb` and
`create-blastdb-file` that work on collections of fasta sequences (see
clj-fasta) and fasta formatted files respectively.

## License

Copyright © 2016 Jason Mulvenna

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
