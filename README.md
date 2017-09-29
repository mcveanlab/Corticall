CortexJDK: A Java API for manipulating (Mc)Cortex de novo assembly graph and link data.
=========

Quick start
-----------

    git clone https://github.com/mcveanlab/CortexJDK
    cd CortexJDK
    ant
    java -jar dist/cortexjdk.jar


Introduction
------------

CortexJDK is a Java class library for performing efficient, low-memory traversals on multi-color linked de Bruijn graphs (LdBG) produced by McCortex.  The most important functionalities provided are:

* iterating over records in a Cortex graphs
* random access (by binary search) to Cortex graph records
* performing simple walks (i.e. extracting a contig, optionally using links to disambiguate junction choices)
* performing walks assisted by one or more reference sequences in a manner consistent with link information
* performing depth-first searches with custom stopping rules (useful for finding interesting graph motifs)
* aligning kmers and contigs back to reference sequences

CortexJDK handles the heavy lifting when operating with these data structures, permitting developers to concentrate on the genome analysis and variant calling tools that can be written on top of this API.

Additionally CortexJDK also contains the graph-based de novo mutation calling software developed during my D.Phil.  These tools may be split into their own repository at a later date.  For now, they serve as a useful example for other developers on how to use the API.


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

To get help for a specific command (e.g. "Print"):

	java -jar dist/cortexjdk.jar Print


Support
-------

Please contact Kiran Garimella (<kiran@well.ox.ac.uk>) with any questions/comments/concerns/cake.  Feedback, bug reports, and pull requests are welcome.


Citing CortexJDK
----------------

A manuscript is currently being prepared.  Other Cortex-related publications are as follows:

* Integrating long-range connectivity information into de Bruijn graphs, Turner, Garimella, Iqbal, McVean (Biorxiv Preprint) (2017) (doi:10.1101/147777) http://www.biorxiv.org/content/early/2017/06/09/147777

* De novo assembly and genotyping of variants using colored de Bruijn graphs, Iqbal, Caccamo, Turner, Flicek, McVean (Nature Genetics) (2012) (doi:10.1038/ng.1028) http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3272472
