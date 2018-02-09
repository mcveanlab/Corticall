package uk.ac.ox.well.cortexjdk.utils.containers;

import com.google.common.base.Joiner;
import uk.ac.ox.well.cortexjdk.utils.io.utils.LineReader;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;

public class DataFrame<R, C, D> {
    private LinkedHashMap<R, LinkedHashMap<C, D>> data = new LinkedHashMap<>();

    private LinkedHashSet<C> colNames = new LinkedHashSet<>();

    private D zeroValue;

    public DataFrame(D zeroValue) {
        this.zeroValue = zeroValue;
    }

    public DataFrame(File d, D zeroValue) {
        this.zeroValue = zeroValue;

        LineReader lr = new LineReader(d);
        String[] header = lr.getNextRecord().split("\t");

        String line;
        while ((line = lr.getNextRecord()) != null) {
            String[] fields = line.split("\t");

            R rowName = (R) fields[0];

            for (int i = 1; i < fields.length; i++) {
                C colName = (C) header[i];
                String entry = fields[i];
                D value = zeroValue;

                if (zeroValue instanceof Float) {
                    value = (D) Float.valueOf((String) entry);
                }

                set(rowName, colName, value);
            }
        }
    }

    public void set(R rowName, C colName, D datum) {
        colNames.add(colName);

        if (!data.containsKey(rowName)) {
            data.put(rowName, new LinkedHashMap<>());
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

    public void addMatrix(DataFrame<R, C, D> d) {
        Collection<R> rows = d.getRowNames();
        Collection<C> cols = d.getColNames();

        for (R r : rows) {
            for (C c : cols) {
                D v0 = this.get(r, c);
                D v1 = d.get(r, c);

                if (v0 instanceof Float && v1 instanceof Float) {
                    Float v2 = (Float) v0 + (Float) v1;
                    this.set(r, c, (D) v2);
                }
            }
        }
    }

    public String toString() {
        String result = "\t" + Joiner.on("\t").join(colNames) + "\n";

        for (R rowName : data.keySet()) {
            ArrayList<String> rowFields = new ArrayList<>();
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
