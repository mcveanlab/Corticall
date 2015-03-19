package uk.ac.ox.well.indiana.utils.containers;

import com.google.common.base.Joiner;

import java.io.PrintStream;
import java.util.*;

public class DataTable implements Iterable<Map<String, Object>>, Iterator<Map<String, Object>> {
    private String tableName;
    private String description;
    private Set<String> columnNames = new LinkedHashSet<String>();
    private Map<String, Map<String, Object>> data = new LinkedHashMap<String, Map<String, Object>>();

    private Iterator<String> pkIterator;

    public DataTable(String tableName, String description) {
        this.tableName = tableName;
        this.description = description;
    }

    public DataTable(String tableName, String description, String... columnNames) {
        this.tableName = tableName;
        this.description = description;

        addColumns(columnNames);
    }

    public void addColumns(String... columnNames) {
        for (String columnName : columnNames) {
            addColumn(columnName);
        }
    }

    public void addColumn(String columnName) {
        columnNames.add(columnName);
    }

    public String getTableName() { return tableName; }
    public String getDescription() { return description; }
    public Set<String> getPrimaryKeys() { return data.keySet(); }

    public boolean has(String primaryKey) {
        return data.containsKey(primaryKey);
    }

    public boolean has(String primaryKey, String columnName) {
        return has(primaryKey) && data.get(primaryKey).containsKey(columnName);
    }

    public void set(String primaryKey, String columnName, Object value) {
        if (!data.containsKey(primaryKey)) {
            data.put(primaryKey, new HashMap<String, Object>());
        }

        if (!columnNames.contains(columnName)) {
            columnNames.add(columnName);
        }

        data.get(primaryKey).put(columnName, value);
    }

    public Object get(String primaryKey, String columnName) {
        //return has(primaryKey, columnName) ? data.get(primaryKey).get(columnName) : "NA";
        return data.get(primaryKey).get(columnName);
    }

    public void write(PrintStream out) {
        out.println(tableName + "\t" + Joiner.on("\t").join(columnNames));

        for (String primaryKey : data.keySet()) {
            List<String> fields = new ArrayList<String>();
            fields.add(tableName);

            for (String columnName : columnNames) {
                fields.add(String.valueOf(get(primaryKey, columnName)));
            }

            out.println(Joiner.on("\t").join(fields));
        }

        out.println();
    }

    @Override
    public Iterator<Map<String, Object>> iterator() {
        pkIterator = data.keySet().iterator();

        return this;
    }

    @Override
    public boolean hasNext() {
        return pkIterator.hasNext();
    }

    @Override
    public Map<String, Object> next() {
        String pk = pkIterator.next();

        return data.get(pk);
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }
}
