package uk.ac.ox.well.indiana.utils.io.utils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.HashMap;
import java.util.Iterator;

public class TableReader2 implements Iterator<HashMap<String, String>>, Iterable<HashMap<String, String>> {
    public TableReader2(File fileToRead) {
        /*
        FileInputStream f = null;
        try {
            f = new FileInputStream(fileToRead);

            FileChannel ch = f.getChannel();
            MappedByteBuffer mb = ch.map( FileChannel.MapMode.READ_ONLY, 0L, ch.size() );

            byte[] barray = new byte[SIZE];
            long checkSum = 0L;
            int nGet;
            while( mb.hasRemaining() ) {
                nGet = Math.min( mb.remaining( ), SIZE );
                mb.get( barray, 0, nGet );
                for ( int i=0; i<nGet; i++ )
                    checkSum += barray[i];
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        */
    }

    @Override
    public Iterator<HashMap<String, String>> iterator() {
        return null;
    }

    @Override
    public boolean hasNext() {
        return false;
    }

    @Override
    public HashMap<String, String> next() {
        return null;
    }

    @Override
    public void remove() {
    }
}
