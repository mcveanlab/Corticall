package uk.ac.ox.well.indiana.utils.arguments;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.cli.*;
import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlEngine;
import uk.ac.ox.well.indiana.commands.Command;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.collection.CortexCollection;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexGraphLinks;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.xmfa.XMFASequenceFile;

import java.awt.*;
import java.io.*;
import java.lang.annotation.Annotation;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.*;
import java.util.List;

public class ArgumentHandler {
    private static JexlEngine je;

    private ArgumentHandler() {}

    public static void parse(Command instance, String[] args) {
        try {
            Options options = new Options();

            Field[] instanceFields = instance.getClass().getDeclaredFields();
            Field[] superFields = instance.getClass().getSuperclass().getDeclaredFields();

            List<Field> fields = new ArrayList<>();
            fields.addAll(Arrays.asList(instanceFields));
            fields.addAll(Arrays.asList(superFields));

            int numArgFields = 0;

            for (Field field : fields) {
                for (Annotation annotation : field.getDeclaredAnnotations()) {
                    if (annotation.annotationType().equals(Argument.class)) {
                        numArgFields++;

                        Argument arg = (Argument) annotation;

                        ArrayList<String> descElements = new ArrayList<>();

                        if (field.get(instance) != null) {
                            descElements.add("default: " + field.get(instance));
                        }

                        if (Collection.class.isAssignableFrom(field.getType()) || Map.class.isAssignableFrom(field.getType())) {
                            descElements.add("can be specified more than once");
                        }

                        if (!arg.required()) {
                            descElements.add("optional");
                        }

                        String description = arg.doc();
                        if (descElements.size() > 0) {
                            description += " [" + Joiner.on(", ").join(descElements) + "]";
                        }

                        Boolean hasArgument = !field.getType().equals(Boolean.class);

                        Option option = new Option(arg.shortName(), arg.fullName(), hasArgument, description);
                        option.setType(field.getType());
                        options.addOption(option);
                    } else if (annotation.annotationType().equals(Output.class)) {
                        numArgFields++;

                        // Todo: this redundancy is fugly

                        Output out = (Output) annotation;

                        String description = out.doc() + " [default: " + (field.getType().equals(PrintStream.class) ? "/dev/stdout" : "/dev/null") + "]";
                        Boolean hasArgument = !field.getType().equals(Boolean.class);

                        Option option = new Option(out.shortName(), out.fullName(), hasArgument, description);
                        option.setType(field.getType());
                        options.addOption(option);
                    }
                }
            }

            options.addOption("h", "help", false, "Show this help message");

            CommandLineParser parser = new IndianaParser();
            CommandLine cmd = parser.parse(options, args);

            if (cmd.hasOption("help") || (numArgFields > 0 && args.length == 0)) {
                int width = System.getenv("COLUMNS") == null ? 100 : Integer.valueOf(System.getenv("COLUMNS"));

                Description d = instance.getClass().getAnnotation(Description.class);
                String header = (d == null) ? "no description available" : d.text();
                String command = "java -jar indiana.jar " + instance.getClass().getSimpleName() + " [arguments]";
                String footer = "";

                HelpFormatter formatter = new HelpFormatter();
                formatter.setWidth(width);
                formatter.setSyntaxPrefix("Usage: ");

                System.out.println();
                formatter.printHelp(command, header, options, footer, false);
                System.out.println();

                System.exit(1);
            }

            for (Field field : fields) {
                for (Annotation annotation : field.getDeclaredAnnotations()) {
                    if (annotation.annotationType().equals(Argument.class)) {
                        Argument arg = (Argument) annotation;

                        if (field.getType().equals(Boolean.class)) {
                            if (cmd.hasOption(arg.fullName())) {
                                Boolean prevValue = (Boolean) field.get(instance);
                                field.set(instance, !prevValue);
                            }
                        } else {
                            String value = cmd.getOptionValue(arg.fullName());

                            if (value != null) {
                                processArgument(instance, field, value);
                            } else if (arg.required() && field.get(instance) == null) {
                                throw new RuntimeException("The argument '--" + arg.fullName() + "' was not specified and is required");
                            }
                        }
                    } else if (annotation.annotationType().equals(Output.class)) {
                        Output out = (Output) annotation;

                        if (field.getType().equals(PrintStream.class) || field.getType().equals(File.class)) {
                            String value = cmd.getOptionValue(out.fullName());

                            if (value == null) {
                                if (field.getType().equals(PrintStream.class)) {
                                    value = "/dev/stdout";
                                } else {
                                    value = "/dev/null";
                                }
                            }

                            processArgument(instance, field, value);
                        } else {
                            String value = cmd.getOptionValue(out.fullName());

                            processArgument(instance, field, value);
                        }
                    }
                }
            }
        } catch (ParseException e) {
            throw new RuntimeException("Error when parsing command-line arguments: " + e);
        } catch (IllegalAccessException e) {
            throw new RuntimeException("Error when accessing field: " + e);
        }
    }

    private static void processArgument(Command instance, Field field, String value) {
        try {
            Class<?> type = field.getType();

            if (Collection.class.isAssignableFrom(type)) {
                ArrayList<String> values = new ArrayList<>();

                for (String avalue : value.split(",")) {
                    File valueAsFile = new File(avalue);
                    if (valueAsFile.exists() && (valueAsFile.getAbsolutePath().endsWith(".list"))) {
                        BufferedReader reader = new BufferedReader(new FileReader(valueAsFile));

                        String line;
                        while ((line = reader.readLine()) != null) {
                            values.add(line);
                        }
                    } else {
                        values.add(avalue);
                    }
                }

                Object o = field.getType().newInstance();

                String containerType = field.getGenericType().toString();
                String genericType = containerType.substring(containerType.indexOf("<") + 1, containerType.lastIndexOf(">"));

                for (String v : values) {
                    Method add = Collection.class.getDeclaredMethod("add", Object.class);
                    add.invoke(o, handleArgumentTypes(Class.forName(genericType), v));
                }

                field.set(instance, o);
            } else if (Map.class.isAssignableFrom(type) && !CortexLinks.class.isAssignableFrom(type)) {
                Map<String, String> pairs = new HashMap<>();

                for (String avalue : value.split(",")) {
                    File valueAsFile = new File(avalue);
                    if (valueAsFile.exists() && (valueAsFile.getAbsolutePath().endsWith(".list") || valueAsFile.getAbsolutePath().endsWith(".txt"))) {
                        BufferedReader reader = new BufferedReader(new FileReader(valueAsFile));

                        String line;
                        while ((line = reader.readLine()) != null) {
                            String[] keyvalue = line.split("[:\\s]+");
                            pairs.put(keyvalue[0], keyvalue[1]);
                        }
                    } else {
                        String[] keyvalue = avalue.split("[:\t]");

                        if (keyvalue.length == 2) {
                            pairs.put(keyvalue[0], keyvalue[1]);
                        } else {
                            throw new IndianaException("Argument for map types must be formatted as key value pairs");
                        }
                    }
                }

                Object o = field.getType().newInstance();

                String containerType = field.getGenericType().toString();
                String genericTypes = containerType.substring(containerType.indexOf("<") + 1, containerType.lastIndexOf(">"));
                String[] keyvalueTypes = genericTypes.replaceAll("\\s+", "").split(",");

                for (String k : pairs.keySet()) {
                    String v = pairs.get(k);

                    Method put = Map.class.getDeclaredMethod("put", Object.class, Object.class);
                    //put.invoke(o, k, handleArgumentTypes(Class.forName(keyvalueTypes[1]), v));
                    put.invoke(o, handleArgumentTypes(Class.forName(keyvalueTypes[0]), k), handleArgumentTypes(Class.forName(keyvalueTypes[1]), v));
                }

                field.set(instance, o);
            } else {
                field.set(instance, handleArgumentTypes(type, value));
            }
        } catch (InvocationTargetException | ClassNotFoundException | NoSuchMethodException | InstantiationException | IllegalAccessException | IOException e) {
            throw new RuntimeException(e);
        }
    }

    private static Object handleArgumentTypes(Class<?> type, String value) {
        try {
            if (type.equals(Integer.class)) {
                return Integer.valueOf(value);
            } else if (type.equals(Long.class)) {
                return Long.valueOf(value);
            } else if (type.equals(Float.class)) {
                return Float.valueOf(value);
            } else if (type.equals(Double.class)) {
                return Double.valueOf(value);
            } else if (type.equals(String.class)) {
                return value;
            } else if (type.equals(CortexKmer.class)) {
                return new CortexKmer(value);
            } else if (type.equals(CortexGraph.class)) {
                return new CortexGraph(value);
            } else if (type.equals(CortexGraphLinks.class)) {
                return new CortexGraphLinks(value);
            } else if (type.equals(CortexLinks.class)) {
                return new CortexLinks(value);
            } else if (type.equals(GFF3.class)) {
                return new GFF3(value);
            } else if (type.equals(FastaSequenceFile.class)) {
                return new FastaSequenceFile(new File(value), false);
            } else if (type.equals(XMFASequenceFile.class)) {
                return new XMFASequenceFile(value);
            } else if (type.equals(IndexedFastaSequenceFile.class)) {
                return new IndexedFastaSequenceFile(new File(value));
            } else if (type.equals(KmerLookup.class)) {
                return new KmerLookup(new File(value));
            } else if (type.equals(PrintStream.class)) {
                FileOutputStream fdout = new FileOutputStream(value);
                BufferedOutputStream bos = new BufferedOutputStream(fdout, 1048576);
                return new PrintStream(bos, false);
            } else if (type.equals(SAMFileReader.class)) {
                SAMFileReader sfr = new SAMFileReader(new File(value));
                sfr.setValidationStringency(ValidationStringency.SILENT);

                return sfr;
            } else if (type.equals(SamReader.class)) {
                return SamReaderFactory.make()
                        .setOption(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, true)
                        .setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, true)
                        .open(new File(value));
            } else if (type.equals(FastqReader.class)) {
                return new FastqReader(new File(value));
            } else if (type.equals(VCFFileReader.class)) {
                return new VCFFileReader(new File(value), false);
            } else if (type.equals(Color.class)) {
                return value.startsWith("#") ? Color.decode(value) : Color.decode("#" + value);
            } else if (type.equals(CortexGraphWriter.class)) {
                return new CortexGraphWriter(value);
            } else if (type.equals(CortexCollection.class)) {
                return new CortexCollection(value);
            } else if (type.equals(Expression.class)) {
                initializeJexlEngine();

                return je.createExpression(value);
            } else if (type.equals(File.class)) {
                return new File(value);
            } else {
                throw new RuntimeException("Unable to automatically handle argument of type '" + type.getSimpleName() + "'");
            }
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Unable to process argument of type '" + type.getSimpleName() + "' and value '" + value + "'", e);
        }
    }

    private static void initializeJexlEngine() {
        if (je == null) {
            je = new JexlEngine();
            je.setCache(512);
            je.setLenient(false);
            je.setSilent(false);
        }
    }
}
