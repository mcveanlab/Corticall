package uk.ac.ox.well.cortexjdk.utils.io.graph.links;

import uk.ac.ox.well.cortexjdk.CortexJDK;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.ConnectivityAnnotations;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexHeader;

import java.io.File;

public class CortexLinks implements ConnectivityAnnotations {
    ConnectivityAnnotations links;

    public CortexLinks(String linksPath) { initialize(new File(linksPath)); }

    public CortexLinks(File linksFile) { initialize(linksFile); }

    private void initialize(File linksFile) {
        File linksIndex = new File(linksFile.getAbsolutePath() + ".idx");

        if (linksIndex.exists()) {
            links = new CortexLinksRandomAccess(linksFile);
        } else {
            //links = new CortexLinksMap(linksFile);

            throw new CortexJDKException("Loading of unindexed links files ('" + linksFile.getAbsolutePath() + "') is disabled");
        }
    }

    @Override
    public int size() {
        return links.size();
    }

    @Override
    public boolean isEmpty() {
        return links.isEmpty();
    }

    @Override
    public boolean containsKey(Object key) {
        return links.containsKey(key);
    }

    @Override
    public CortexLinksRecord get(Object key) {
        return links.get(key);
    }

    @Override
    public CortexHeader getHeader() {
        return links.getHeader();
    }

    @Override
    public String getSource() { return links.getSource(); }
}
