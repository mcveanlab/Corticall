package uk.ac.ox.well.indiana.utils.packageutils;

import uk.ac.ox.well.indiana.tools.Tool;

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
            Class<? extends Tool> module = (Class<? extends Tool>) Class.forName(moduleName);

            Tool instance = module.newInstance();
            instance.args = moduleArgs;

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
