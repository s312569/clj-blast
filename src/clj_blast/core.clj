(ns clj-blast.core
  (:require [fs.core :refer [absolute-path temp-file delete]]
            [clojure.data.xml :refer [parse]]
            [clojure.zip :refer [xml-zip node]]
            [clojure.data.zip.xml :refer [text xml1-> xml->]]
            [clj-commons-exec :refer [sh]]
            [clj-fasta.core :refer [fasta->file fasta-seq]]
            [clojure.java.io :refer [reader]]))

(defn iteration-seq
  "Returns a lazy list of iterations in a blast XML file."
  [reader]
  (->> (parse reader)
       :content
       (filter #(= (:tag %) :BlastOutput_iterations))
       first
       :content
       (filter #(= (:tag %) :Iteration))
       (map xml-zip)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; accessors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-hsp-value
  "Takes a blastHsp object and returns the value corresponding to key.
   Keys are the keyword version of the XML nodes in the BLAST xml
   output.  All values are returned as strings. Typical BLAST HSP
   keys are:
    :Hsp_bit-score
    :Hsp_score
    :Hsp_evalue
    :Hsp_query-from
    :Hsp_query-to
    :Hsp_hit-from
    :Hsp_hit-to
    :Hsp_positive
    :Hsp_identity
    :Hsp_gaps
    :Hsp_hitgaps
    :Hsp_querygaps
    :Hsp_qseq
    :Hsp_hseq
    :Hsp_midline
    :Hsp_align-len
    :Hsp_query-frame
    :Hsp_hit-frame
    :Hsp_num
    :Hsp_pattern-from
    :Hsp_pattern-to
    :Hsp_density"
  [hsp key]
  (and hsp
       (if (= key :Hsp_midline)
         (->> (:content (node hsp))
              (filter #(= (:tag %) :Hsp_midline))
              first :content first)
         (xml1-> hsp key text))))

(defn get-hit-value
  "Takes a blastHit object and returns the value corresponding to key.
  Keys are the keyword version of the XML nodes in the BLAST xml
  output.  All values are returned as strings. Typical BLAST Hit
  values are:
   :Hit_id
   :Hit_len
   :Hit_accession
   :Hit_def
   :Hit_num"
  [hit key]
  (xml1-> hit key text))

(defn hit-seq
  "Returns a lazy list of hits (as zippers)."
  [it]
  (and it (xml-> it :Iteration_hits :Hit)))

(defn hsp-seq
  "Returns a lazy list of HSPs (as zippers)."
  [hit]
  (and hit (xml-> hit :Hit_hsps :Hsp)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; alignment
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

(defn hsp-alignment
  "Returns the alignment from a HSP as a string."
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
          (throw (Throwable.
                  (str "Blast error: " (:err bl))))
          (throw (Throwable.
                  (str "Exception: " (:exception bl)))))))))

(defn blast
  "Takes a collection of fasta sequences and blasts them against the
  specified database."
  [bs program db outfile & {:keys [params] :or {params {}}}]
  (let [i (absolute-path (fasta->file bs (temp-file "seq-")))]
    (try
      (run-blast program (absolute-path db) i (absolute-path outfile) params)
      (finally (delete i)))))

(defn blast-file
  "Blasts a file with fasta formatted sequences. Blasts 10,000
  sequences at a time in parallel using pmap."
  [file program db outfile & {:keys [params] :or {params {}}}]
  (let [c (atom 0)]
    (with-open [r (reader file)]
      (doall
       (pmap #(blast % program db (str outfile "-" (swap! c inc) ".xml"))
             (partition-all 10000 (fasta-seq r)))))))
