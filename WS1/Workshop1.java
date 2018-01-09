/*
 * Compile: javac *.java
 * Run: java Workshop1UI
 */
public class Workshop1 {
	/**
	 * Task 1:
	 * Retrieve the intensity value at location (x, y) of the image and return it
	 * Note: 
	 * - the 2D image is stored as an 8bit, 1D, row-major array of type byte
	 * - Note that byte is a signed data type in Java
	 * @param img in row major format
	 * @param x coordinate
	 * @param y coordinate
	 * @param width of img
	 * @param height of img
	 * @return the intensity value at (x, y) if within bounds, -1 otherwise
	 */
	public int getIntensityValue(byte[] img, int x, int y, int width, int height) {
		System.out.println("w ="+width);
		System.out.println("h ="+height);
		return (int)(img[x*width+y] & 0xFF);
	}
	/**
	 * Task 2:
	 * Retrieve the intensity value that occurs most often in the image
	 * @param img 
	 * @return the intensity value that occurs most often in the image
	 */
	public int getMostFrequentIntensityValue(byte[] img) {
		int[] f = new int[256];
		for (int ib : img)
			f[(int)(ib & 0xFF)]++;
		int max = -1;
		for (int i=0;i<256;i++)
			max  = (max < f[i]) ? i : max;
		return max;
	}	

	private void changeABunchOfPixels(byte[] img, int width, int[] coordinates, int newValue) {
		for (int i=0;i<coordinates.length;i+=2)
			img[coordinates[i]*width+coordinates[i+1]] = (byte)newValue;
	}

	public void setEightNeighborsToWhite(byte[] img, int x, int y, int width, int height) {
		//System.out.println("TODO: implement Task 3");
		int [] pixelsToBeChanged = new int[]{x-1, y-1, x-1, y, x-1, y+1, 
			x, y-1, x, y+1, x+1, y-1, x+1, y, x+1, y+1};
		changeABunchOfPixels(img, width, pixelsToBeChanged, 255);

		/*int [] pixelsToBeChanged2 = new int[100*50*2];
		int i=0;
		for (int xc = 150;xc<250;xc++) {
			for (int yc = 250;yc<300;yc++) {
				pixelsToBeChanged2[i++] = xc;
				pixelsToBeChanged2[i++] = yc;
			}
		}
		changeABunchOfPixels(img, width, pixelsToBeChanged2, 120);*/
	}
	/**
	 * Task 4:
	 * Calculates the d4 distance between (x1, y1) and (x2, y2)
	 * @param img that will be unchanged
	 * @param x1 
	 * @param y1
	 * @param x2
	 * @param y2
	 * @param width of img
	 * @param height of img
	 * @return the d4 distance between (x1, y1) and (x2, y2) 
	 */
	public int getD4Distance(byte[] img, int x1, int y1, int x2, int y2, int width, int height) {
		return Math.abs(x1-x2)+Math.abs(y1-y2);
	}
	/**
	 * Homework 1: 
	 * Marks the shortest m-path with white intensity values. Let V = {0, ..., 127}.
	 * This task was developed to challenge yourself. Can you find a shortest m-path quickly?
	 * @param img with the shortest m-path set to white with V = {0, ..., 127}.
	 * @param x1 
	 * @param y1
	 * @param x2
	 * @param y2
	 * @param width of img
	 * @param height of img
	 */
	public void setShortestMPathToWhite(byte[] img, int x1, int y1, int x2, int y2, int width, int height) {
		System.out.println("TODO: implement Homework 1");
	}
}
