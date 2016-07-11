(defproject clj-blast "0.2.4"
  :description "Parser for BLAST output in XML format."
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.8.0"]
                 [org.clojure/tools.nrepl "0.2.12"]
                 [org.clojure/data.xml "0.0.8"]
                 [org.clojure/data.zip "0.1.1"]
                 [clj-fasta "0.1.7"]
                 [biodb "0.1.9"]
                 [org.clojars.hozumi/clj-commons-exec "1.2.0"]
                 [fs "1.3.3"]]
  :repl-options {:init (set! *print-length* 100)}
  :jvm-opts ["-Xmx1000M"])
