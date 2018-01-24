package uk.ac.ox.well.cortexjdk.utils.traversal;

public class LinkStoreElement {
    private String junctionList;
    private int age;
    private int pos;
    private String source;

    public LinkStoreElement(String junctionList, int age, int pos, String source) {
        this.junctionList = junctionList;
        this.age = age;
        this.pos = pos;
        this.source = source;
    }

    public void incrementAge() { age++; }

    public void incrementPos() { pos++; }

    public boolean hasAgedOut() { return age >= junctionList.length(); }

    public String getJunctionList() { return junctionList; }

    public int getAge() { return age; }

    public int getPos() { return pos; }

    public String getSource() { return source; }

    public int length() { return junctionList.length(); }

    @Override
    public String toString() {
        return junctionList + " [" + pos + "/" + junctionList.length() + "] age: " + age;
    }
}
