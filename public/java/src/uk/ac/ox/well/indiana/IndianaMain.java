package uk.ac.ox.well.indiana;

import org.apache.commons.cli.*;
import uk.ac.ox.well.indiana.exceptions.IndianaException;
import uk.ac.ox.well.indiana.exceptions.code.ToolResolutionException;
import uk.ac.ox.well.indiana.exceptions.user.ToolNotFoundException;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.packageutils.PackageInspector;

import java.lang.annotation.Annotation;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;

public class IndianaMain {
    public static void main(String[] args) throws IndianaException {
        if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help") || args[0].equals("-H")) {
            showHelp();
        } else if (args.length > 0) {
            String toolName = args[0];
            String[] toolArgs = Arrays.copyOfRange(args, 1, args.length);

            HashMap<String, Class<? extends Tool>> tools = new PackageInspector<Tool>(Tool.class).getExtendingClassesMap();

            if (tools.containsKey(toolName)) {
                Class<? extends Tool> tool = tools.get(toolName);

                try {
                    Tool instance = tool.newInstance();
                    Options options = new Options();

                    Field[] fields = instance.getClass().getDeclaredFields();

                    for (Field field : fields) {
                        for (Annotation annotation : field.getDeclaredAnnotations()) {
                            if (annotation.annotationType().equals(Argument.class)) {
                                Argument arg = (Argument) annotation;

                                options.addOption(arg.shortName(), arg.fullName(), arg.required(), arg.doc());
                            }
                        }
                    }

                    CommandLineParser parser = new PosixParser();
                    CommandLine cmd = parser.parse(options, toolArgs);

                    /*
                    HelpFormatter formatter = new HelpFormatter();
                    formatter.printHelp("java -jar indiana.jar " + toolName, options);

                    Option[] parsedOptions = cmd.getOptions();
                    for (Option option : parsedOptions) {
                        System.out.println(option + " " + option.getValue());
                    }
                    */

                    for (Field field : fields) {
                        for (Annotation annotation : field.getDeclaredAnnotations()) {
                            if (annotation.annotationType().equals(Argument.class)) {
                                Argument arg = (Argument) annotation;

                                String value = cmd.getOptionValue(arg.fullName());
                                field.set(instance, value);
                            }
                        }
                    }

                    Method method = instance.getClass().getMethod("execute");

                    method.invoke(instance);
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
                throw new ToolNotFoundException(toolName);
            }
        }
    }

    private static void showHelp() {
        Set<Class<? extends Tool>> tools = new PackageInspector<Tool>(Tool.class).getExtendingClasses();

        System.out.println("usage: java -jar indiana.jar [-h|--help|-H]");
        System.out.println("                             <command> [<args>]");
        System.out.println();

        System.out.println("tools:");
        for (Class t : tools) {
            System.out.println("   " + t.getSimpleName());
        }
        System.out.println();
    }
}