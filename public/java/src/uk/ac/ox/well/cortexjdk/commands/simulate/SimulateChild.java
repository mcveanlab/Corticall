package uk.ac.ox.well.cortexjdk.commands.simulate;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;

/**
 * Created by kiran on 19/11/2017.
 */
public class SimulateChild extends Module {
    @Argument(fullName="ref1", shortName="r1", doc="Ref 1")
    public IndexedFastaSequenceFile REF1;

    @Argument(fullName="ref2", shortName="r2", doc="Ref 2")
    public IndexedFastaSequenceFile REF2;

    @Override
    public void execute() {

}
