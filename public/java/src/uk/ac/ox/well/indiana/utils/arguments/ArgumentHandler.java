package uk.ac.ox.well.indiana.utils.arguments;

import com.google.common.base.Joiner;
import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.apache.commons.cli.*;
import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlEngine;
import uk.ac.ox.well.indiana.IndianaModule;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;

import java.io.*;
import java.lang.annotation.Annotation;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

public class ArgumentHandler {
    private static JexlEngine je;

    private ArgumentHandler() {}

    public static void parse(IndianaModule instance, String[] args) {
        try {
            Options options = new Options();

            Field[] instanceFields = instance.getClass().getDeclaredFields();
            Field[] superFields = instance.getClass().getSuperclass().getDeclaredFields();

            List<Field> fields = new ArrayList<Field>();
            fields.addAll(Arrays.asList(instanceFields));
            fields.addAll(Arrays.asList(superFields));

            int numArgFields = 0;

            for (Field field : fields) {
                for (Annotation annotation : field.getDeclaredAnnotations()) {
                    if (annotation.annotationType().equals(Argument.class)) {
                        numArgFields++;

                        Argument arg = (Argument) annotation;

                        ArrayList<String> descElements = new ArrayList<String>();

                        if (field.get(instance) != null) {
                            descElements.add("default: " + field.get(instance));
                        }

                        if (Collection.class.isAssignableFrom(field.getType())) {
                            descElements.add("can be specified more than once");
                        }

                        String description = arg.doc();
                        if (descElements.size() > 0) {
                            description += " [" + Joiner.on(" ").join(descElements) + "]";
                        }

                        Boolean hasArgument = !field.getType().equals(Boolean.class);

                        Option option = new Option(arg.shortName(), arg.fullName(), hasArgument, description);
                        option.setType(field.getType());
                        options.addOption(option);
                    } else if (annotation.annotationType().equals(Output.class)) {
                        numArgFields++;

                        // Todo: this redundancy is fugly

                        Output out = (Output) annotation;

                        String description = out.doc() + " [default: /dev/stdout]";
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
                HelpFormatter formatter = new HelpFormatter();

                int width = System.getenv("COLUMNS") == null ? 100 : Integer.valueOf(System.getenv("COLUMNS"));

                formatter.setWidth(width);
                formatter.printHelp("java -jar indiana.jar " + instance.getClass().getSimpleName() + " [options]", options);
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

                        if (field.getType().equals(PrintStream.class)) {
                            String value = cmd.getOptionValue(out.fullName());

                            if (value == null) {
                                value = "/dev/stdout";
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

    private static void processArgument(IndianaModule instance, Field field, String value) {
        try {
            Class<?> type = field.getType();

            if (Collection.class.isAssignableFrom(type)) {
                ArrayList<String> values = new ArrayList<String>();

                for (String avalue : value.split(",")) {
                    File valueAsFile = new File(avalue);
                    if (valueAsFile.exists() && valueAsFile.getAbsolutePath().endsWith(".list")) {
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
            } else {
                field.set(instance, handleArgumentTypes(type, value));
            }
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        } catch (IOException e) {
            throw new RuntimeException(e);
        } catch (IllegalAccessException e) {
            throw new RuntimeException(e);
        } catch (InstantiationException e) {
            throw new RuntimeException(e);
        } catch (NoSuchMethodException e) {
            throw new RuntimeException(e);
        } catch (ClassNotFoundException e) {
            throw new RuntimeException(e);
        } catch (InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    private static Object handleArgumentTypes(Class<?> type, String value) {
        try {
            if (type.equals(File.class)) {
                return new File(value);
            } else if (type.equals(Integer.class)) {
                return Integer.valueOf(value);
            } else if (type.equals(Float.class)) {
                return Float.valueOf(value);
            } else if (type.equals(Double.class)) {
                return Double.valueOf(value);
            } else if (type.equals(String.class)) {
                return value;
            } else if (type.equals(CortexGraph.class)) {
                return new CortexGraph(value);
            } else if (type.equals(GFF3.class)) {
                return new GFF3(value);
            } else if (type.equals(FastaSequenceFile.class)) {
                return new FastaSequenceFile(new File(value), false);
            } else if (type.equals(IndexedFastaSequenceFile.class)) {
                return new IndexedFastaSequenceFile(new File(value));
            } else if (type.equals(PrintStream.class)) {
                return new PrintStream(value);
            } else if (type.equals(Expression.class)) {
                initializeJexlEngine();

                return je.createExpression(value);
            } else {
                throw new RuntimeException("Unable to automatically handle argument of type '" + type.getSimpleName() + "'");
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        return null;
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
