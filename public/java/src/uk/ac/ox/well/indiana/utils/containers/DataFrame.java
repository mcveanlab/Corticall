package uk.ac.ox.well.indiana.utils.containers;

import com.google.common.base.Joiner;

import java.util.*;

public class DataFrame<R, C, D> {
    private LinkedHashMap<R, LinkedHashMap<C, D>> data = new LinkedHashMap<R, LinkedHashMap<C, D>>();

    private LinkedHashSet<C> colNames = new LinkedHashSet<C>();

    private D zeroValue;

    public DataFrame(D zeroValue) {
        this.zeroValue = zeroValue;
    }

    public void set(R rowName, C colName, D datum) {
        colNames.add(colName);

        if (!data.containsKey(rowName)) {
            data.put(rowName, new LinkedHashMap<C, D>());
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

    public boolean hasValue(R rowName, C colName) {
        return (data.containsKey(rowName) && data.get(rowName).containsKey(colName));
    }

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
