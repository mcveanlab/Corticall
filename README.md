CortexJDK
=========

Author: Kiran Garimella <kiran@well.ox.ac.uk>

A library for manipulating (Mc)Cortex de-novo assembly graph and link data.

Dependencies
------------

You must have Apache Ant (http://ant.apache.org) installed to build CortexJDK. Other dependencies are fetched via Ivy at compile time.

McCortex is required for building assemblies and creating link annotations.

Installation
------------

To download dependencies and compile

    git clone https://github.com/kvg/CortexJDK.git
    cd CortexJDK
    ant

To run:

    java -jar dist/cortexjdk.jar
