package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class ClipOverlappingContigs extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file of aligned contigs")
    public SAMFileReader BAM;

    @Output
    public File out;

    private boolean isDovetail(SAMRecord lastContig, SAMRecord thisContig) {
        return (lastContig.getAlignmentEnd() > thisContig.getAlignmentStart() && lastContig.getAlignmentEnd() < thisContig.getAlignmentEnd());
    }

    private boolean isContained(SAMRecord lastContig, SAMRecord thisContig) {
        return (thisContig.getAlignmentStart() >= lastContig.getAlignmentStart() && thisContig.getAlignmentEnd() <= lastContig.getAlignmentEnd());
    }

    private int computeOverlapAmount(SAMRecord lastContig, SAMRecord thisContig) {
        //return lastContig.getReadLength() - (thisContig.getAlignmentStart() - lastContig.getAlignmentStart());
        return lastContig.getAlignmentEnd() - thisContig.getAlignmentStart() + 1;
    }

    private Cigar getNonRedundantCigar(List<CigarElement> ces) {
        List<CigarElement> newCes = new ArrayList<CigarElement>();
        CigarOperator curOp = ces.get(0).getOperator();
        int curLength = ces.get(0).getLength();

        for (int i = 1; i < ces.size(); i++) {
            if (curOp.equals(ces.get(i).getOperator())) {
                curLength += ces.get(i).getLength();
            } else {
                newCes.add(new CigarElement(curLength, curOp));

                curOp = ces.get(i).getOperator();
                curLength = ces.get(i).getLength();
            }
        }

        newCes.add(new CigarElement(curLength, curOp));

        return new Cigar(newCes);
    }

    private String getCigarString(Cigar cigar) {
        StringBuilder cigarString = new StringBuilder();

        for (CigarElement ce : cigar.getCigarElements()) {
            cigarString.append(ce.getLength()).append(ce.getOperator().toString());
        }

        return cigarString.toString();
    }

    private int numOps(List<CigarElement> ces, CigarOperator op) {
        int numOps = 0;

        for (CigarElement ce : ces) {
            if (op.equals(ce.getOperator())) {
                numOps++;
            }
        }

        return numOps;
    }

    @Override
    public void execute() {
        SAMFileWriter sfw = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(BAM.getFileHeader(), false, out);

        log.info("Processing contigs...");
        int numAlignmentsChanged = 0, numAlignmentsRemoved = 0, numAlignments = 0;
        SAMRecord lastContig = null;

        for (SAMRecord thisContig : BAM) {
            if (thisContig.getMappingQuality() > 0) {
                if (lastContig != null && lastContig.getReferenceName().equals(thisContig.getReferenceName())) {
                    if (isDovetail(lastContig, thisContig)) {
                        int overlapAmount = computeOverlapAmount(lastContig, thisContig);

                        //log.debug("{}", lastContig.getReadString());
                        //log.debug("{}{}", StringUtil.repeatCharNTimes(' ', thisContig.getAlignmentStart() - lastContig.getAlignmentStart()), thisContig.getReadString());
                        //log.debug("{}", overlapAmount);

                        int overlapUsed = 0;
                        int alignmentShift = 0;
                        List<CigarElement> oldCigarList = thisContig.getCigar().getCigarElements();
                        List<CigarElement> newCigarList = new ArrayList<CigarElement>();
                        for (int i = 0; i < oldCigarList.size(); i++) {
                            CigarElement ce = oldCigarList.get(i);
                            CigarOperator co = ce.getOperator();

                            if (overlapUsed < overlapAmount) {
                                if (co.equals(CigarOperator.H) || co.equals(CigarOperator.S)) {
                                    newCigarList.add(ce);
                                } else if (co.equals(CigarOperator.M)) {
                                    if (ce.getLength() <= overlapAmount - overlapUsed) {
                                        newCigarList.add(new CigarElement(ce.getLength(), CigarOperator.S));

                                        overlapUsed += ce.getLength();
                                    } else {
                                        newCigarList.add(new CigarElement(overlapAmount - overlapUsed, CigarOperator.S));
                                        newCigarList.add(new CigarElement(ce.getLength() - (overlapAmount - overlapUsed), CigarOperator.M));

                                        overlapUsed = overlapAmount;
                                    }
                                } else if (co.equals(CigarOperator.I)) {
                                    newCigarList.add(new CigarElement(ce.getLength(), CigarOperator.S));
                                    overlapUsed += ce.getLength();

                                    alignmentShift -= ce.getLength();
                                } else if (co.equals(CigarOperator.D)) {
                                    // do nothing?
                                    alignmentShift += ce.getLength();
                                }
                            } else {
                                newCigarList.add(ce);
                            }
                        }

                        //Cigar cigar = new Cigar(newCigarList);
                        Cigar newCigar = getNonRedundantCigar(newCigarList);

                        thisContig.setCigar(newCigar);
                        thisContig.setAlignmentStart(thisContig.getAlignmentStart() + overlapAmount + alignmentShift);

                        //log.debug("{} {} {} {}", lastContig.getCigarString(), thisContig.getCigarString(), cigar, getCigarString(newCigar));
                        //log.debug("-----");

                        numAlignmentsChanged++;
                    } else if (isContained(lastContig, thisContig)) {
                        thisContig.setMappingQuality(0);
                        thisContig.setReadUnmappedFlag(true);

                        numAlignmentsRemoved++;
                    }
                }

                sfw.addAlignment(thisContig);

                lastContig = thisContig;
            } else {
                sfw.addAlignment(thisContig);
            }

            numAlignments++;
        }

        sfw.close();

        log.info("  changed: {}/{}", numAlignmentsChanged, numAlignments);
        log.info("  removed: {}/{}", numAlignmentsRemoved, numAlignments);
    }
}
