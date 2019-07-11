CortexJDK: A Java API for manipulating (Mc)Cortex de novo assembly graph and link data.
=========

[![Build Status](https://travis-ci.org/mcveanlab/CortexJDK.svg?branch=master)](https://travis-ci.org/mcveanlab/CortexJDK)

Quick start
-----------

    git clone https://github.com/mcveanlab/CortexJDK
    cd CortexJDK
    ant
    java -jar dist/cortexjdk.jar -h


Introduction
------------

CortexJDK is a Java class library for performing efficient, low-memory traversals on multi-color linked de Bruijn graphs (LdBG) produced by McCortex.  The most important functionalities provided are:

* iterating over records in a Cortex graphs
* random access (by binary search) to Cortex graph records
* performing simple walks (i.e. extracting a contig, optionally using links to disambiguate junction choices)
* performing walks assisted by one or more reference sequences in a manner consistent with link information
* performing depth-first searches with custom stopping rules (useful for finding interesting graph motifs)
* aligning k-mers and contigs back to reference sequences

CortexJDK handles the heavy lifting when operating with these data structures, permitting developers to concentrate on the genome analysis and variant calling tools that can be written on top of this API.

CortexJDK also contains the graph-based de novo mutation calling software, Corticall.  This tool may be split into its own repository at a later date.


Availability
------------

CortexJDK is released under the [Apache 2.0](https://opensource.org/licenses/Apache-2.0) license.  The latest code is [freely available at Github](https://github.com/mcveanlab/CortexJDK).


Dependencies
------------

CortexJDK has the following dependencies:

* [Java8](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html): needed for runtime and development kit
* [Apache Ant](http://ant.apache.org): for dependency fetching and compilation
* [McCortex](https://github.com/mcveanlab/mccortex): for building Cortex graphs and link annotations

Other dependencies are automatically fetched at compile time.

Installation
------------

To download and compile:

    git clone https://github.com/mcveanlab/CortexJDK
    cd CortexJDK
    ant

To run (and get a listing of available commands):

    java -jar dist/cortexjdk.jar

To get help for a specific command (e.g. "Call"):

    java -jar dist/cortexjdk.jar Call


Support
-------

Please contact Kiran V Garimella (<kiran@broadinstitute.org>) with any questions/comments/concerns/cake.  Feedback, bug reports, and pull requests are welcome.


Related projects
----------------
For a pythonic take on programmatically accessing McCortex linked multi-color de Bruijn graphs, see Winni Kretzschmar's [cortexpy](https://pypi.org/project/cortexpy/) project.


Citation
--------

The original multi-color de Bruijn graph paper can be found at:

* "De novo assembly and genotyping of variants using colored de Bruijn graphs", Iqbal, Caccamo, Turner, Flicek, McVean (Nature Genetics) (2012) (doi:10.1038/ng.1028) http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3272472

Our manuscript on linked multi-color de Bruijn graphs is now available via Bioinformatics at:

* "Integrating long-range connectivity information into de Bruijn graphs", Turner, Garimella, Iqbal, McVean (Bioinformatics) (2018) (doi:10.1093/bioinformatics/bty157) https://www.ncbi.nlm.nih.gov/pubmed/29554215

A manuscript preprint for Corticall, our de novo assembly approach to de novo mutation detection in pathogens, can be found at:

* "Detection of simple and complex de novo mutations without, with, or with multiple reference sequences", Garimella, Iqbal, Krause, Campino, Kekre, Drury, Kwiatkowski, Sa, Wellems, McVean (2019) (doi:10.1101/698910) https://www.biorxiv.org/content/10.1101/698910v1

A manuscript for CortexJDK is currently being prepared.
