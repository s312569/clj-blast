# clj-blast

A parser for BLAST XML files.

## Usage

Import from Clojars:

```clojure
[clj-blast "0.1.2"]
```

Use in your namespace:

```clojure
(require '[clj-blast.core :as bl])
```

Open a reader on a blast output file and call 'iteration-seq' on the
reader to get a lazy list of iterations. To get hits call 'hit-seq' on
iterations and likewise 'hsp-seq' on hits to get hsps. Hits and hsps
have 'get-hit-value' and 'get-hsp-value' for accessing data. 

```clojure
user> (with-open [r (reader "/blast/output/file.xml")]
        (get-hsp-value (->> (iteration-seq r)
                            first
                            hit-seq
                            first
                            hsp-seq
                            first)
                       :Hsp_bit-score))
"23.8682"
user>
```

Blasts can be performed using 'blast' and 'blast-file'.

## License

Copyright © 2016 Jason Mulvenna

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
