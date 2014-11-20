/**
 * Copyright 2009-2010 by The Broad Institute.
 *
 */
package cmap.matrix;

/**
 * Default metadata implementation.
 */
public class DefaultMetadataModel extends AbstractMetadataModel {

    private int itemCount;

    public DefaultMetadataModel(int itemCount) {
	this.itemCount = itemCount;
    }

    @Override
    public int getItemCount() {
	return itemCount;
    }

    @Override
    protected Vector createVector(String name, Class<?> c) {
	return new DefaultVector(name, c, getItemCount());
    }

}