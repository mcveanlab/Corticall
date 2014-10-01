INDIANA
=======

Author: Kiran Garimella <kiran@well.ox.ac.uk>

Tools to manipulate and analyze genome assemblies.

Dependencies
------------

You must have Apache Ant (http://ant.apache.org) installed to build INDIANA. Other dependencies are fetched via Ivy at compile time.

Installation
------------

To download dependencies and compile

    git clone https://github.com/mcveanlab/INDIANA.git
    cd INDIANA
    ant

To run:

    java -jar dist/indiana.jar

Usage
-----

View a cortex graph:

    ctx31 sort graph.ctx
    samtools faidx seq.fa
    java -jar dist/indiana.jar ctxgraph -p seq.k9.ctp.gz -ieg seq.k9.ctx -c seq.fa -o seq.dot
    dot -Tpdf seq.dot > seq.pdf

LICENSE
-------

GPLv2
