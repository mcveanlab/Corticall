package uk.ac.ox.well.indiana.utils.arguments;

import net.sf.picard.reference.FastaSequenceFile;
import org.apache.commons.cli.*;
import uk.ac.ox.well.indiana.IndianaModule;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;

import java.io.*;
import java.lang.annotation.Annotation;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

public class ArgumentParser {
    private ArgumentParser() {}

    public static void parse(IndianaModule instance, String[] args) {
        try {
            Options options = new Options();

            Field[] fields = instance.getClass().getDeclaredFields();

            for (Field field : fields) {
                for (Annotation annotation : field.getDeclaredAnnotations()) {
                    if (annotation.annotationType().equals(Argument.class)) {
                        Argument arg = (Argument) annotation;

                        String description = arg.doc() + " [default: " + field.get(instance) + "]";
                        Boolean hasArgument = !field.getType().equals(Boolean.class);

                        Option option = new Option(arg.shortName(), arg.fullName(), hasArgument, description);
                        option.setType(field.getType());
                        options.addOption(option);
                    } else if (annotation.annotationType().equals(Output.class)) {
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

            CommandLineParser parser = new PosixParser();
            CommandLine cmd = parser.parse(options, args);

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
                                handleArgumentTypes(instance, field, value);
                            } else if (arg.required()) {
                                // do something else here
                            }
                        }
                    } else if (annotation.annotationType().equals(Output.class)) {
                        Output out = (Output) annotation;

                        if (field.getType().equals(PrintStream.class)) {
                            String value = cmd.getOptionValue(out.fullName());

                            if (value == null) {
                                value = "/dev/stdout";
                            }

                            handleArgumentTypes(instance, field, value);
                        }
                    }
                }
            }

            if (cmd.hasOption("help") || args.length == 0) {
                HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp("java -jar indiana.jar " + instance.getClass().getSimpleName() + " [options]", options);
                System.exit(1);
            }
        } catch (ParseException e) {
            throw new RuntimeException("Error when parsing command-line arguments: " + e);
        } catch (IllegalAccessException e) {
            throw new RuntimeException("Error when accessing field: " + e);
        }
    }

    /*
    private static void handleArgumentTypes(IndianaModule instance, Field field, String value) {
       try {
           Class<?> type = field.getType();

           //if (type.isArray()) {
           if (Collection.class.isAssignableFrom(type)) {
               ArrayList<String> values = new ArrayList<String>();

               File valueAsFile = new File(value);
               if (valueAsFile.exists() && valueAsFile.getAbsolutePath().endsWith(".list")) {
                   BufferedReader reader = new BufferedReader(new FileReader(valueAsFile));

                   String line;
                   while ((line = reader.readLine()) != null) {
                       values.add(line);
                   }
               } else if (value.contains(",")) {
                   values.addAll(Arrays.asList(value.split(",")));
               }



               System.out.println(type.getComponentType());
               System.out.println(arrayInstance.getClass().getComponentType());

               for (int i = 0; i < arrayInstance.length; i++) {
                   //arrayInstance[i] = (Object) values.get(i);
               }

               field.set(instance, arrayInstance);
           } else {

           }
       } catch (FileNotFoundException e) {
           e.printStackTrace();
       } catch (IOException e) {
           e.printStackTrace();
       } catch (IllegalAccessException e) {
           e.printStackTrace();
       }
    }
    */

    private static void handleArgumentTypes(IndianaModule instance, Field field, String value) {
        try {
            Class<?> type = field.getType();

            if (type.equals(File.class)) {
                field.set(instance, new File(value));
            } else if (type.equals(Integer.class)) {
                field.set(instance, Integer.valueOf(value));
            } else if (type.equals(Float.class)) {
                field.set(instance, Float.valueOf(value));
            } else if (type.equals(Double.class)) {
                field.set(instance, Double.valueOf(value));
            } else if (type.equals(String.class)) {
                field.set(instance, value);
            } else if (type.equals(CortexGraph.class)) {
                field.set(instance, new CortexGraph(value));
            } else if (type.equals(FastaSequenceFile.class)) {
                field.set(instance, new FastaSequenceFile(new File(value), false));
            } else if (type.equals(PrintStream.class)) {
                field.set(instance, new PrintStream(value));
            } else {
                throw new RuntimeException("Don't know how to automatically handle field type '" + type.getSimpleName() + "' for field '" + field.getName() + "'");
            }
        } catch (IllegalAccessException e) {
            throw new RuntimeException("Unable to access field '" + field.getName() + "' in module '" + instance.getClass().getSimpleName() + "'");
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Unable to open file: " + e);
        }
    }
}
