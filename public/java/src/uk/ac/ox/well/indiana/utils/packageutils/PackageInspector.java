package uk.ac.ox.well.indiana.utils.packageutils;

import org.reflections.Reflections;
import org.reflections.scanners.SubTypesScanner;
import org.reflections.util.ClasspathHelper;
import org.reflections.util.ConfigurationBuilder;

import java.lang.reflect.Modifier;
import java.net.URL;
import java.util.*;

/**
 * Inspects the package by reflection and finds classes that extend a specific class.
 *
 * @param <ClassType>  the type of class to look for in the package
 */
public class PackageInspector<ClassType> {
    private Reflections reflections;
    private Class<ClassType> classType;
    private String packageName;

    /**
     * Create a package inspector that looks for the specified type.
     *
     * @param classType  the type of class to look for in the package
     */
    public PackageInspector(Class<ClassType> classType, String packageName) {
        Collection<URL> classPath = ClasspathHelper.forPackage(packageName);

        ConfigurationBuilder config = new ConfigurationBuilder();
        config.setUrls(classPath);
        config.setScanners(new SubTypesScanner());

        this.reflections = new Reflections(config);
        this.classType = classType;
        this.packageName = packageName;
    }

    /**
     * Get the set of classes that extend the specified type.
     *
     * @return  the set of classes that extend the type
     */
    public Set<Class<? extends ClassType>> getExtendingClasses() {
        Set<Class<? extends ClassType>> extendingClasses = reflections.getSubTypesOf(classType);
        Set<Class<? extends ClassType>> nonAbstractClasses = new HashSet<>();

        for (Class<? extends ClassType> c : extendingClasses) {
            if (!Modifier.isAbstract(c.getModifiers())) {
                if (c.getPackage().getName().contains(packageName)) {
                    nonAbstractClasses.add(c);
                }
            }
        }

        return nonAbstractClasses;
    }

    /**
     * Get the map of class name to class object that extend the specified type.
     *
     * @return  the set of classes that extend the type
     */
    public Map<String, Class<? extends ClassType>> getExtendingClassesMap() {
        Map<String, Class<? extends ClassType>> classHashMap = new TreeMap<>();

        for (Class<? extends ClassType> c : getExtendingClasses()) {
            classHashMap.put(c.getSimpleName(), c);
        }

        return classHashMap;
    }

    /**
     * Get the tree of package and class name to class object that extend the specified type.
     *
     * @return  the tree of classes that extend the type
     */
    public Map<String, Map<String, Class<? extends ClassType>>> getExtendingClassTree() {
        Map<String, Map<String, Class<? extends ClassType>>> classTree = new TreeMap<>();

        for (Class<? extends ClassType> c : getExtendingClasses()) {
            String baseName = c.getPackage().getName().replaceAll(packageName + ".", "");
            String className = c.getSimpleName();

            if (!classTree.containsKey(baseName)) {
                classTree.put(baseName, new TreeMap<>());
            }

            classTree.get(baseName).put(className, c);
        }

        return classTree;
    }
}
