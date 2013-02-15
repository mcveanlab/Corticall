package uk.ac.ox.well.indiana.utils.containers;

import com.google.common.base.Joiner;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.TreeSet;

public class DataFrame<R, C, D> {
    private HashMap<R, HashMap<C, D>> data = new HashMap<R, HashMap<C, D>>();

    private TreeSet<C> colNames = new TreeSet<C>();

    private D zeroValue;

    public DataFrame(D zeroValue) {
        this.zeroValue = zeroValue;
    }

    public void set(R rowName, C colName, D datum) {
        colNames.add(colName);

        if (!data.containsKey(rowName)) {
            data.put(rowName, new HashMap<C, D>());
        }

        if (!data.get(rowName).containsKey(colName)) {
            data.get(rowName).put(colName, zeroValue);
        }

        data.get(rowName).put(colName, datum);
    }

    public Collection<R> getRowNames() {
        return data.keySet();
    }

    public Collection<C> getColNames() {
        return colNames;
    }

    public int getNumRows() { return data.size(); }
    public int getNumCols() { return colNames.size(); }

    public D get(R rowName, C colName) {
        if (data.containsKey(rowName) && data.get(rowName).containsKey(colName)) {
            return data.get(rowName).get(colName);
        }

        return zeroValue;
    }

    public String toString() {
        String result = "\t" + Joiner.on("\t").join(colNames) + "\n";

        for (R rowName : data.keySet()) {
            ArrayList<String> rowFields = new ArrayList<String>();
            rowFields.add(rowName.toString());

            for (C colName : colNames) {
                D datum = get(rowName, colName);

                rowFields.add(datum.toString());
            }

            result += Joiner.on("\t").join(rowFields) + "\n";
        }

        return result;
    }
}
