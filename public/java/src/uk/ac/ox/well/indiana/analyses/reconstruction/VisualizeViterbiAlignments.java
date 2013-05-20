package uk.ac.ox.well.indiana.analyses.reconstruction;

import uk.ac.ox.well.indiana.sketches.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;

import java.io.File;

public class VisualizeViterbiAlignments extends Sketch {
    @Argument(fullName="alignment", shortName="a", doc="Alignment file")
    public File ALIGNMENT;

    @Argument(fullName="kmerReferencePanel", shortName="krp", doc="Kmer reference panel")
    public File KMER_REFERENCE_PANEL;

    public void initialize() {
    }
}
