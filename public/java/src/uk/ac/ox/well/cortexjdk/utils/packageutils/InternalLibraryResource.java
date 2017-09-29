package uk.ac.ox.well.cortexjdk.utils.packageutils;

import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.Field;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.Arrays;

/**
 * Created by kiran on 29/09/2017.
 */
public class InternalLibraryResource {
    private final File tempResource;

    public InternalLibraryResource(String resourcePath) {
        File resourceFile = new File(resourcePath);

        try {
            File tempDir = Files.createTempDirectory("bwalib").toFile();

            tempResource = new File(tempDir.getAbsolutePath() + "/" + resourceFile.getName());
            tempResource.deleteOnExit();

            if (!tempResource.exists()) {
                InputStream is = this.getClass().getResourceAsStream(resourcePath);
                Files.copy(is, tempResource.toPath(), StandardCopyOption.REPLACE_EXISTING);
            }

            addLibraryPath(tempDir.getAbsolutePath());
        } catch (IOException e) {
            throw new CortexJDKException("Could not copy resource to temp dir");
        } catch (Exception e) {
            throw new CortexJDKException("Could not set custom library path");
        }
    }

    public File getFile() { return tempResource; }

    // From: http://fahdshariff.blogspot.co.uk/2011/08/changing-java-library-path-at-runtime.html
    public static void addLibraryPath(String pathToAdd) throws Exception{
        final Field usrPathsField = ClassLoader.class.getDeclaredField("usr_paths");
        usrPathsField.setAccessible(true);

        //get array of paths
        final String[] paths = (String[])usrPathsField.get(null);

        //check if the path to add is already present
        for(String path : paths) {
            if(path.equals(pathToAdd)) {
                return;
            }
        }

        //add the new path
        final String[] newPaths = Arrays.copyOf(paths, paths.length + 1);
        newPaths[newPaths.length-1] = pathToAdd;
        usrPathsField.set(null, newPaths);
    }
}
