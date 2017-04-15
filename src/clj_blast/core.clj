(ns clj-blast.core
  (:require [me.raynes.fs :as fs]
            [clojure.data.xml :refer [parse]]
            [clojure.zip :refer [xml-zip node]]
            [clojure.data.zip.xml :refer [text xml1-> xml->]]
            [clj-commons-exec :refer [sh]]
            [biodb.core :as bdb]
            [clojure.string :refer [split trim]]
            [clj-fasta.core :refer [fasta->file fasta-seq]]
            [clojure.java.io :refer [reader]]))

(declare hit->hm hsp-alignment)

(defn- get-hit-value
  [hit key]
  (xml1-> hit key text))

(defn- get-hsp-value
  [hsp key]
  (and hsp
       (if (= key :Hsp_midline)
         (->> (:content (node hsp))
              (filter #(= (:tag %) :Hsp_midline))
              first :content first)
         (xml1-> hsp key text))))

(defn- hsp->hm
  [hsp]
  (let [doubles [:Hsp_bit-score :Hsp_score :Hsp_evalue]
        integers [:Hsp_query-from :Hsp_query-to :Hsp_hit-from :Hsp_hit-to 
                  :Hsp_pattern-to :Hsp_query-frame :Hsp_hit-frame
                  :Hsp_identity :Hsp_positive :Hsp_gaps :Hsp_align-len
                  :Hsp_density :Hsp_pattern-from :Hsp_num]
        others [:Hsp_qseq :Hsp_hseq :Hsp_midline]]
    (-> (merge (into {} (map #(vector % (if-let [v (get-hsp-value hsp %)] (Double/parseDouble v)))
                             doubles))
               (into {} (map #(vector % (if-let [v (get-hsp-value hsp %)] (Integer/parseInt v)))
                             integers))
               (into {} (map #(vector % (get-hsp-value hsp %)) others)))
        (assoc :alignment (hsp-alignment hsp)))))

(defn- hit->hm
  [hit]
  (-> (->> (map #(vector % (get-hit-value hit %))
                [:Hit_id :Hit_len :Hit_accession :Hit_def :Hit_num])
           (into {}))
      (assoc :hsps (map hsp->hm (xml-> hit :Hit_hsps :Hsp)))
      (update-in [:Hit_len] #(if % (Integer/parseInt %)))
      (update-in [:Hit_num] #(if % (Integer/parseInt %)))))

(defn- residue-counter [start end]
  (let [c (atom start)]
    {:increment (if (< end start)
                  (fn [x] (reset! c (- @c x)))
                  (fn [x] (reset! c (+ @c x))))
     :value (fn [] (deref c))}))

(defn- formatter
  [start end]
  (let [c (residue-counter start end)]
    (fn [x]
      (let [n (count (remove #{\-} x))]
        (vector (str ((:value c)))
                (apply str x)
                (str (- ((:increment c) n) 1)))))))

(defn- get-lines
  [hsp keys]
  (let [f (partial get-hsp-value hsp)
        args (map f keys)
        form (apply formatter (map #(Integer/parseInt %)
                                   (drop 1 args)))]
    (map form (partition-all 52 (first args)))))

(defn- hsp-alignment
  [hsp]
  (let [l (interleave
           (get-lines hsp '(:Hsp_qseq :Hsp_query-from :Hsp_query-to))
           (map #(vector "" (apply str %) "")
                (partition-all 52 (get-hsp-value hsp :Hsp_midline)))
           (get-lines hsp '(:Hsp_hseq :Hsp_hit-from :Hsp_hit-to)))
        m (apply max (mapcat #(list (count (first %))
                                    (count (nth % 2)))
                             l))
        b (fn [s] (apply str
                         (repeat (+ 2 (- m (count s))) \space)))
        pf (fn [x] (str (first x) (b (first x)) (second x) "  "
                        (last x) "\n"))]
    (apply str (interpose "\n" (map #(apply str %)
                                    (partition-all 3 (map pf l)))))))

(defn- significant?
  [hit param score]
  (condp = param
    :evalue (some #(<= (:Hsp_evalue %) score) (:hsps hit))
    :bit-score (some #(>= (:Hsp_bit-score %) score) (:hsps hit))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; accessors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn iteration-seq
  "Returns a lazy list of zippers representing iterations in a blast
  XML file."
  [reader]
  (->> (parse reader)
       :content
       (filter #(= (:tag %) :BlastOutput_iterations))
       first
       :content
       (filter #(= (:tag %) :Iteration))
       (map xml-zip)))

(defn accession
  "Returns the accession of the query sequence in an iteration."
  [it]
  (and it (-> (xml1-> it :Iteration_query-def text)
              (split #"\s+")
              first)))

(defn query-length
  "Returns the query length from an iteration."
  [it]
  (if-let [ql (and it (xml1-> it :Iteration_query-len text))]
    (Integer/parseInt ql)))

(defn query-def
  "Returns the def line of the query from an iteration."
  [it]
  (and it (xml1-> it :Iteration_query-def text)))

(defn iteration-number
  "Returns the iteration number."
  [it]
  (if-let [in (and it (xml1-> it :Iteration_iter-num text))]
    (Integer/parseInt in)))

(defn hit-seq
  "Returns a lazy list of hits as maps representing a blast hit. HSPs
  can be accessed using :hsps. Alignments can be accessed
  using :alignment on individual HSP maps. If an evalue or bit-score
  keyword argument is supplied hits will be filtered for
  significance."
  [it & {:keys [evalue bit-score] :or [evalue nil bit-score nil]}]
  (if it
    (let [hits (map #(assoc (hit->hm %) :query-accession (accession it))
                    (xml-> it :Iteration_hits :Hit))]
      (cond (and evalue bit-score) (throw (Exception. "Supply one of :bit-score or :evalue only."))
            evalue (filter #(significant? % :evalue evalue) hits)
            bit-score (filter #(significant? % :bit-score bit-score) hits)
            :else hits))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; blasting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- blast-default-params
  [params in-file out-file db]
  (let [p (dissoc params "-ungapped")
        ungapped? (if (get params "-ungapped")
                    "-ungapped")]
    (doall
     (->> (merge {"-evalue" "10"
                  "-outfmt" "5"
                  "-max_target_seqs" "3"
                  "-query"
                  in-file
                  "-out"
                  out-file
                  "-db" db}
                 p)
          seq
          flatten
          (cons ungapped?)
          (remove nil?)))))

(defn- run-blast 
  [prog db in out params]
  (let [defs (blast-default-params params
                                   in
                                   out
                                   db)]
    (let [bl @(sh (cons prog defs))]
      (if (= 0 (:exit bl))
        out
        (if (:err bl)
          (throw (Exception.
                  (str "Blast error: " (:err bl))))
          (throw (Exception.
                  (str "Exception: " (:exception bl)))))))))

(defn- blast-partition
  [bs program db outfile & {:keys [params] :or {params {}}}]
  (let [i (fs/absolute (fasta->file bs (fs/temp-file "seq-")))]
    (try
      (run-blast program (fs/absolute db) i (fs/absolute outfile) params)
      (finally (fs/delete i)))))

(defn blast
  "Takes a collection of fasta sequences and blasts them against the
  specified database. Blasts 10,000 sequences at a time in parallel
  using pmap."
  [bs program db outfile & {:keys [params] :or {params {}}}]
  (let [c (atom 0)]
    (doall
     (->> (pmap #(blast-partition % program db
                                  (str (fs/base-name outfile true) "-" (swap! c inc)
                                       (fs/extension outfile))
                                  :params params)
                (partition-all 10000 bs))
          flatten))))

(defn blast-file
  "Blasts a file with fasta formatted sequences. Blasts 10,000
  sequences at a time in parallel using pmap."
  [file program db outfile & {:keys [params] :or {params {}}}]
  (let [c (atom 0)]
    (with-open [r (reader file)]
      (blast (fasta-seq r) program db outfile :params params))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; sequence retrieval
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- blastdbcommand-cmd
  ([coll db dbtype] (blastdbcommand-cmd coll db dbtype nil))
  ([coll db dbtype outfile]
   (let [accs (if (> (count coll) 1000)
                ["-entry_batch" (str (fasta->file coll (fs/temp-file "bdbc-") :func identity))]
                ["-entry" (->> (map trim coll) (interpose ",") (apply str))])
         out (if outfile ["-out" (str outfile)] [])]
     (try
       (let [rbs @(sh (concat ["blastdbcmd" "-db" (str db) "-dbtype" dbtype] accs out))]
         (if (= (:exit rbs) 0)
           (or outfile
               (with-open [r (java.io.BufferedReader. (java.io.StringReader. (:out rbs)))]
                 (doall (fasta-seq r))))
           (if (:err rbs)
             (if-not (re-find #"OID not found" (:err rbs))
               (throw (Exception. (str "Blast error: " (:err rbs)))))
             (throw (Exception. (str "Exception: " (:exception rbs)))))))
       (finally
         (if (string? accs)
           (fs/delete accs)))))))

(defn retrieve-sequence
  "Takes a collection of accession numbers and retrieves a collection
  of fasta sequences from a blast database. NOT LAZY!"
  ([coll db] (retrieve-sequence coll db "prot"))
  ([coll db dbtype]
   (blastdbcommand-cmd coll db dbtype)))

(defn blastdb->file
  "Takes a collection of accessions, retrieves them from a blastdb and
  sends them to the outfile in fasta format."
  ([coll db outfile] (blastdb->file coll db outfile "prot"))
  ([coll db outfile dbtype]
   (blastdbcommand-cmd coll db dbtype outfile)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; db creation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn create-blastdb
  "Takes a collection of fasta sequences, writes them to the specified
  file and creates a blast database. Returns the blast database path."
  ([coll db-name] (create-blastdb coll db-name "prot"))
  ([coll db-name dbtype]
   (let [faf (fasta->file coll db-name)]
     (try
       (let [db @(sh ["makeblastdb" "-in" (str faf) "-dbtype" dbtype "-parse_seqids"])]
         (if (= (:exit db) 0)
           faf
           (if (:err db)
             (throw (Exception. (str "makeblastdb error: " (:err db))))
             (throw (Exception. (str "Exception: " (:exception db)))))))))))

(defn create-blastdb-file
  "Takes a path to a file of fasta sequences and creates a blast
  database. Returns the blast database path."
  ([file] (create-blastdb-file file "prot"))
  ([file dbtype]
   (try
     (let [db @(sh ["makeblastdb" "-in" (str file) "-dbtype" dbtype "-parse_seqids"])]
       (if (= (:exit db) 0)
         file
         (if (:err db)
           (throw (Exception. (str "makeblastdb error: " (:err db))))
           (throw (Exception. (str "Exception: " (:exception db))))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; biodb compatability
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod bdb/table-spec :blast
  [q]
  (vector [:accession :text "PRIMARY KEY"]
          [:src :binary "NOT NULL"]))

(defmethod bdb/prep-sequences :blast
  [q]
  (->> (:coll q)
       (map #(hash-map :accession (accession %) :src (bdb/freeze (node %))))))

(defmethod bdb/restore-sequence :blast
  [q]
  (xml-zip (bdb/thaw (:src (dissoc q :type)))))

