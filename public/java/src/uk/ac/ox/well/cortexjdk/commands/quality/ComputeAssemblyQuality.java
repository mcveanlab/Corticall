package uk.ac.ox.well.cortexjdk.commands.quality;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.collection.CortexCollection;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableWriter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by kiran on 08/08/2017.
 */
public class ComputeAssemblyQuality extends Module {
    @Argument(fullName="eval", shortName="e", doc="Graph to evaluate")
    public CortexGraph EVAL;

    @Argument(fullName="comp", shortName="c", doc="Graph to compare against (i.e. truth)")
    public CortexGraph COMP;

    @Argument(fullName="evalRef", shortName="r", doc="")
    public KmerLookup REF;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        CortexCollection cc = new CortexCollection(EVAL, COMP);

        Set<CortexKmer> seeds = getVariantSeeds(cc, 0, 1);

        int numBases = 0;
        for (SAMSequenceRecord ssr : REF.getReferenceSequence().getSequenceDictionary().getSequences()) {
            numBases += ssr.getSequenceLength();
        }

        log.info("Q: {}", -10*Math.log10((double) seeds.size() / (double) numBases));

        out.println(-10*Math.log10((double) seeds.size() / (double) numBases));

        //int numVariants = callVariants(cc, 0, 1, seeds);
        //log.info("Found {} variants", numVariants);
    }

    private int callVariants(CortexCollection cc, int evalColor, int compColor, Set<CortexKmer> seeds) {
        TraversalEngine e = new TraversalEngineFactory()
                .graph(cc)
                .make();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Building contigs")
                .message("seeds processed")
                .maxRecord(seeds.size())
                .make(log);

        Map<Interval, Map<String, String>> entries = new TreeMap<>();

        int numVariants = 0;
        for (CortexKmer ck : seeds) {
            Map<String, String> te = callVariant(cc, evalColor, compColor, e, ck);

            if (te != null) {
                Interval it = new Interval(te.get("chrom"), Integer.valueOf(te.get("pos")), Integer.valueOf(te.get("pos")));
                entries.put(it, te);

                numVariants++;
            }

            pm.update();
        }

        TableWriter tw = new TableWriter(out);

        for (Interval it : entries.keySet()) {
            tw.addEntry(entries.get(it));
        }

        return numVariants;
    }

    private Map<String, String> callVariant(CortexCollection cc, int evalColor, int compColor, TraversalEngine e, CortexKmer ck) {
        CortexRecord cr = cc.findRecord(ck);
        if (cr.getCoverage(evalColor) > 0) {
            e.getConfiguration().setTraversalColor(evalColor);

            String sk = ck.getKmerAsString();
            List<CortexVertex> contigEval = new ArrayList<>();
            contigEval.add(new CortexVertex(new CortexByteKmer(sk), cr));

            CortexVertex source = null;
            e.seek(sk);
            while (e.hasPrevious()) {
                CortexVertex cv = e.previous();
                contigEval.add(0, cv);

                if (cv.getCr().getCoverage(compColor) > 0 && (cv.getSk().equals(cv.getCk().getKmerAsString()) ? cv.getCr().getOutDegree(compColor) == 1 : cv.getCr().getInDegree(compColor) == 1)) {
                    source = cv;
                    break;
                }
            }

            CortexVertex destination = null;
            e.seek(sk);
            while (e.hasNext()) {
                CortexVertex cv = e.next();
                contigEval.add(cv);

                if (cv.getCr().getCoverage(compColor) > 0 && (cv.getSk().equals(cv.getCk().getKmerAsString()) ? cv.getCr().getInDegree(compColor) == 1 : cv.getCr().getOutDegree(compColor) == 1)) {
                    destination = cv;
                    break;
                }
            }

            if (source != null && destination != null) {
                e.getConfiguration().setTraversalColor(compColor);

                List<CortexVertex> contigComp = new ArrayList<>();
                contigComp.add(source);

                boolean destinationReached = false;

                e.seek(source.getSk());
                while (e.hasNext()) {
                    CortexVertex cv = e.next();
                    contigComp.add(cv);

                    if (cv.equals(destination)) {
                        destinationReached = true;
                        break;
                    }
                }

                if (destinationReached) {
                    Set<Interval> sourceIntervals = REF.findKmer(source.getSk());
                    Set<Interval> destinationIntervals = REF.findKmer(destination.getSk());

                    float compCoverage = 0.0f;
                    for (CortexVertex cv : contigComp) {
                        compCoverage += cv.getCr().getCoverage(compColor);
                    }
                    compCoverage /= (float) contigComp.size();

                    if (sourceIntervals.size() == 1 && destinationIntervals.size() == 1) {
                        Interval sourceInterval = sourceIntervals.iterator().next();
                        Interval destinationInterval = destinationIntervals.iterator().next();

                        if (sourceInterval.getContig().equals(destinationInterval.getContig())) {
                            Pair<String, String> alleles = trimToAlleles(TraversalEngine.toContig(contigEval), TraversalEngine.toContig(contigComp));

                            Map<String, String> te = new LinkedHashMap<>();
                            te.put("chrom", sourceInterval.getContig());
                            te.put("pos", String.valueOf(sourceInterval.getStart()));

                            if (alleles.getFirst().length() == 1 && alleles.getSecond().length() == 1) {
                                te.put("type", "SNP");
                            } else if (alleles.getFirst().length() == alleles.getSecond().length()) {
                                te.put("type", "MNP");
                            } else if (alleles.getFirst().length() < alleles.getSecond().length()) {
                                te.put("type", "DEL");
                            } else if (alleles.getFirst().length() > alleles.getSecond().length()) {
                                te.put("type", "INS");
                            } else {
                                te.put("type", "UKN");
                            }

                            float evalCoverage = 0.0f;
                            for (CortexVertex cv : contigEval) {
                                evalCoverage += cv.getCr().getCoverage(evalColor);
                            }
                            evalCoverage /= (float) contigEval.size();

                            te.put("cov_comp", String.valueOf((int) compCoverage));
                            te.put("cov_eval", String.valueOf((int) evalCoverage));
                            te.put("alleles", String.format("%s/%s", alleles.getFirst(), alleles.getSecond()));

                            return te;
                        }
                    }
                }
            }
        }

        return null;
    }

    @NotNull
    private Set<CortexKmer> getVariantSeeds(CortexCollection cc, int evalColor, int compColor) {
        log.info("Finding variant seeds...");

        Set<CortexKmer> seeds = new HashSet<>();

        DirectedGraph<String, DefaultEdge> sg = new DefaultDirectedGraph<>(DefaultEdge.class);

        for (CortexRecord cr : cc) {
            if (isSinglyConnected(cr) &&
                isUniqueToEval(cr, evalColor, compColor) &&
                hasUniqueCoordinates(cr)) {

                seeds.add(cr.getCortexKmer());

                String skFwd = cr.getKmerAsString();
                String skRev = SequenceUtils.reverseComplement(cr.getKmerAsString());
                sg.addVertex(skFwd);
                sg.addVertex(skRev);

                for (int c = 0; c < cr.getNumColors(); c++) {
                    if (cr.getCoverage(c) > 0) {
                        Collection<String> ies = cr.getInEdgesAsStrings(c, false);
                        for (String ie : ies) {
                            String inEdgeFwd = ie + skFwd.substring(0, skFwd.length() - 1);
                            String outEdgeRev = SequenceUtils.reverseComplement(inEdgeFwd);

                            sg.addVertex(inEdgeFwd);
                            sg.addEdge(inEdgeFwd, skFwd);

                            sg.addVertex(outEdgeRev);
                            sg.addEdge(skRev, outEdgeRev);
                        }

                        Collection<String> oes = cr.getOutEdgesAsStrings(c, false);
                        for (String oe : oes) {
                            String outEdgeFwd = skFwd.substring(1, skFwd.length()) + oe;
                            String inEdgeRev = SequenceUtils.reverseComplement(outEdgeFwd);

                            sg.addVertex(outEdgeFwd);
                            sg.addEdge(skFwd, outEdgeFwd);

                            sg.addVertex(inEdgeRev);
                            sg.addEdge(inEdgeRev, skRev);
                        }
                    }
                }
            }
        }

        Set<String> uniqueSeeds = new HashSet<>();
        Set<CortexKmer> goodSeeds = new HashSet<>();
        for (String sk : sg.vertexSet()) {
            if (sg.inDegreeOf(sk) == 0 && sg.outDegreeOf(sk) == 1) {
                uniqueSeeds.add(sk);


                Set<String> contig = new LinkedHashSet<>();
                contig.add(sk);

                String v = sk;
                while (sg.outDegreeOf(v) == 1) {
                    List<String> out = Graphs.successorListOf(sg, v);
                    v = out.get(0);

                    if (contig.contains(v)) {
                        break;
                    }

                    contig.add(v);
                }

                if (contig.size() > 3) {
                    List<String> linContig = new ArrayList<>(contig);
                    goodSeeds.add(new CortexKmer(linContig.get(1)));
                }
            }
        }

        log.info("  found {} seed kmers for putative variants, {} unique, {} good", seeds.size(), uniqueSeeds.size(), goodSeeds.size());

        return goodSeeds;
    }

    private boolean isSinglyConnected(CortexRecord cr) {
        boolean isSinglyConnected = true;
        for (int c = 0; c < cr.getNumColors(); c++) {
            if (cr.getCoverage(c) > 0) {
                if (!(cr.getInDegree(c) == 1 && cr.getOutDegree(c) == 1)) {
                    isSinglyConnected = false;
                }
            }
        }

        return isSinglyConnected;
    }

    private boolean isUniqueToEval(CortexRecord cr, int evalColor, int compColor) {
        return cr.getCoverage(evalColor) > 0 && cr.getCoverage(compColor) == 0;
    }

    private boolean hasUniqueCoordinates(CortexRecord cr) {
        Set<Interval> its = REF.findKmer(cr.getKmerAsString());

        return its != null && its.size() == 1;
    }

    private Pair<String, String> trimToAlleles(String s0, String s1) {
        int s0start = 0, s0end = s0.length();
        int s1start = 0, s1end = s1.length();

        for (int i = 0, j = 0; i < s0.length() && j < s1.length(); i++, j++) {
            if (s0.charAt(i) != s1.charAt(j)) {
                s0start = i;
                s1start = j;
                break;
            }
        }

        for (int i = s0.length() - 1, j = s1.length() - 1; i >= 0 && j >= 0; i--, j--) {
            if (s0.charAt(i) != s1.charAt(j) || i == s0start - 1 || j == s1start - 1) {
                s0end = i + 1;
                s1end = j + 1;
                break;
            }
        }

        String[] pieces = new String[4];
        pieces[0] = s0.substring(0, s0start);
        pieces[1] = s0.substring(s0start, s0end);
        pieces[2] = s1.substring(s1start, s1end);
        pieces[3] = s0.substring(s0end, s0.length());

        return new Pair<>(pieces[1], pieces[2]);
    }
}
