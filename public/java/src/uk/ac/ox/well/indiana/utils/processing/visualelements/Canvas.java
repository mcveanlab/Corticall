package uk.ac.ox.well.indiana.utils.processing.visualelements;

import java.util.ArrayList;
import java.util.List;

public class Canvas implements VisualElement {
    private List<VisualElement> visualElements = new ArrayList<VisualElement>();

    public void addElement(VisualElement e) {
        visualElements.add(e);
    }

    @Override
    public void draw() {
        for (VisualElement e : visualElements) {
            e.draw();
        }
    }

    @Override
    public int getWidth() {
        int width = 0;

        for (VisualElement e : visualElements) {
            if (e.getWidth() > width) {
                width = e.getWidth();
            }
        }

        return width;
    }

    @Override
    public int getHeight() {
        int height = 0;

        for (VisualElement e : visualElements) {
            if (e.getHeight() > height) {
                height = e.getHeight();
            }
        }

        return height;
    }
}
