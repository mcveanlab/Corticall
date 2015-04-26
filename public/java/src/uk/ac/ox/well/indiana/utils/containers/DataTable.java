package uk.ac.ox.well.indiana.utils.containers;

import com.google.common.base.Joiner;

import java.io.PrintStream;
import java.util.*;

public class DataTable implements Iterable<Map<String, Object>>, Iterator<Map<String, Object>> {
    private String tableName;
    private String description;
    private Set<String> columnNames = new LinkedHashSet<String>();
    private Map<String, Map<String, Object>> data = new TreeMap<String, Map<String, Object>>();

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
        return data.get(primaryKey).get(columnName);
    }

    public void increment(String primaryKey, String columnName) {
        Object o = has(primaryKey, columnName) ? get(primaryKey, columnName) : 0l;

        set(primaryKey, columnName, ((Long) o) + 1);
    }

    public void decrement(String primaryKey, String columnName) {
        Object o = has(primaryKey, columnName) ? get(primaryKey, columnName) : 0l;

        set(primaryKey, columnName, ((Long) o) - 1);
    }

    public void add(String primaryKey, String columnName, long o1) {
        Object o = has(primaryKey, columnName) ? get(primaryKey, columnName) : 0l;

        set(primaryKey, columnName, ((Long) o) + o1);
    }

    public void write(PrintStream out) {
        out.println(this.toString());
    }

    @Override
    public String toString() {
        Map<String, Integer> columnWidths = new HashMap<String, Integer>();
        for (String primaryKey : data.keySet()) {
            for (String columnName : columnNames) {
                String field = String.valueOf(get(primaryKey, columnName));

                if (!columnWidths.containsKey(columnName)) {
                    columnWidths.put(columnName, columnName.length());
                }

                int oldLength = columnWidths.get(columnName);

                if (field.length() > oldLength) {
                    columnWidths.put(columnName, field.length());
                }
            }
        }

        StringBuilder sb = new StringBuilder();
        //sb.append(tableName).append("\t").append(Joiner.on("\t").join(columnNames)).append("\n");

        sb.append(tableName);
        for (String columnName : columnNames) {
            int fieldWidth = columnWidths.get(columnName);
            sb.append(" ");
            sb.append(String.format("%-" + fieldWidth + "s", columnName));
        }
        sb.append("\n");

        for (String primaryKey : data.keySet()) {
            List<String> fields = new ArrayList<String>();
            fields.add(tableName);

            for (String columnName : columnNames) {
                String field = String.valueOf(get(primaryKey, columnName));
                int fieldWidth = columnWidths.get(columnName);

                fields.add(String.format("%-" + fieldWidth + "s", field));

                //fields.add(String.valueOf(get(primaryKey, columnName)));
            }

            //sb.append(Joiner.on("\t").join(fields)).append("\n");
            sb.append(Joiner.on(" ").join(fields)).append("\n");
        }

        sb.append("\n");

        return sb.toString();
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
