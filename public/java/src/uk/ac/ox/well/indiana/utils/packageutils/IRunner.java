package uk.ac.ox.well.indiana.utils.packageutils;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;

import java.lang.annotation.Annotation;
import java.lang.reflect.Field;
import java.util.*;

/**
 * Instantiates a module and passes it command-line arguments.
 */
public class IRunner {
    /**
     * Private constructor - this class cannot be instantiated!
     */
    private IRunner() {}

    /**
     * Main method for instantiating INDIANA modules
     *
     * @param moduleName  the name of the module
     * @param moduleArgs  the arguments for the module
     */
    @SuppressWarnings("unchecked")
    public static void main(String moduleName, String[] moduleArgs) {
        try {
            Class<? extends Module> module = (Class<? extends Module>) Class.forName(moduleName);

            Module instance = module.newInstance();
            instance.args = moduleArgs;

            Set<String> specifiedArgs = new HashSet<String>();
            for (int i = 0; i < moduleArgs.length - 1; i+=2) {
                specifiedArgs.add(moduleArgs[i]);
            }

            Field[] instanceFields = instance.getClass().getDeclaredFields();
            Field[] superFields = instance.getClass().getSuperclass().getDeclaredFields();

            List<Field> fields = new ArrayList<Field>();
            fields.addAll(Arrays.asList(instanceFields));
            fields.addAll(Arrays.asList(superFields));

            List<String> defaultArgs = new ArrayList<String>();
            for (Field field : fields) {
                for (Annotation annotation : field.getDeclaredAnnotations()) {
                    if (annotation.annotationType().equals(Argument.class)) {
                        Argument arg = (Argument) annotation;

                        if (!specifiedArgs.contains(arg.shortName()) && !specifiedArgs.contains(arg.fullName())) {
                            if (field.get(instance) != null) {
                                Object defaultValue = field.get(instance);

                                defaultArgs.add("-" + arg.shortName() + " " + defaultValue.toString());
                            }
                        }
                    }
                }
            }

            instance.log.info("Invocation: {} {}    {}", instance.getClass().getSimpleName(), Joiner.on(" ").join(moduleArgs), Joiner.on(" ").join(defaultArgs));
            instance.log.info("");

            instance.init();
            instance.execute();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        } catch (InstantiationException e) {
            e.printStackTrace();
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        }
    }
}
