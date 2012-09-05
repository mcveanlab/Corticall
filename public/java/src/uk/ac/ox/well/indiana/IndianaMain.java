package uk.ac.ox.well.indiana;

import org.apache.commons.cli.*;
import uk.ac.ox.well.indiana.exceptions.IndianaException;
import uk.ac.ox.well.indiana.exceptions.code.ToolResolutionException;
import uk.ac.ox.well.indiana.exceptions.user.ToolNotFoundException;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.packageutils.PackageInspector;

import java.io.File;
import java.lang.annotation.Annotation;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;

public class IndianaMain {
    public static void main(String[] args) throws Exception {
        if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help") || args[0].equals("-H")) {
            showPrimaryHelp();
        } else if (args.length > 0) {
            String toolName = args[0];
            String[] toolArgs = Arrays.copyOfRange(args, 1, args.length);

            HashMap<String, Class<? extends Tool>> tools = new PackageInspector<Tool>(Tool.class).getExtendingClassesMap();

            if (tools.containsKey(toolName)) {
                try {
                    Tool instance = tools.get(toolName).newInstance();
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
                            }
                        }
                    }

                    options.addOption("h", "help", false, "Show this help message");

                    CommandLineParser parser = new PosixParser();
                    CommandLine cmd = parser.parse(options, toolArgs);

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
                                        } else {
                                        }
                                    } else if (arg.required()) {
                                        // do something else here
                                    }
                                }
                            }
                        }
                    }

                    if (cmd.hasOption("help") || toolArgs.length == 0) {
                        HelpFormatter formatter = new HelpFormatter();
                        formatter.printHelp("java -jar indiana.jar " + toolName + " [options]", options);
                        System.exit(1);
                    }

                    Method method = instance.getClass().getMethod("execute");
                    Integer result = (Integer) method.invoke(instance);

                    System.exit(result);
                } catch (InstantiationException e) {
                    throw new ToolResolutionException(toolName, e);
                } catch (IllegalAccessException e) {
                    throw new ToolResolutionException(toolName, e);
                } catch (ParseException e) {
                    e.printStackTrace();
                } catch (NoSuchMethodException e) {
                    e.printStackTrace();
                } catch (InvocationTargetException e) {
                    e.printStackTrace();
                }
            } else {
                showInvalidToolMessage(toolName);
            }
        }
    }

    private static void showPrimaryHelp() {
        Set<Class<? extends Tool>> tools = new PackageInspector<Tool>(Tool.class).getExtendingClasses();

        System.out.println("usage: java -jar indiana.jar [-h|--help|-H]");
        System.out.println("                             <command> [<args>]");
        System.out.println();

        System.out.println("tools:");
        for (Class t : tools) {
            System.out.println("   " + t.getSimpleName());
        }
        System.out.println();

        System.exit(1);
    }

    private static void showInvalidToolMessage(String tool) {
        System.out.println("indiana: '" + tool + "' is not an INDIANA tool. See 'java -jar indiana.jar --help'.");
        System.exit(1);
    }
}