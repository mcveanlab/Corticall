package uk.ac.ox.well.indiana.utils.assembly;

import org.jgrapht.ext.EdgeNameProvider;

public class CortexEdgeNameProvider implements EdgeNameProvider<CortexEdge> {
    @Override
    public String getEdgeName(CortexEdge cortexEdge) {
        return cortexEdge.getLabel();
    }
}
