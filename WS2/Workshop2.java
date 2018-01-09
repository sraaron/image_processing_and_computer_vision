public class Workshop2 {
	/**
	 * Task 1: 
     * Implement the thresholding transformation. It will set 
     * all pixels with intensity below the threshold 
     * to black and other pixels to white.
     * @param img
     * @param threshold
     */
    public void thresholdTransformation(byte[] img, int threshold) {
    	for (int i=0;i<img.length;i++) {
    		if ((int)(img[i] & 0xFF) < threshold)
    			img[i] = (byte)0;
    		else
    			img[i] = (byte)255;
    	}
    }
   /**
     * Task 2: 
	 * Implement the negative transformation.
	 *
	 * @param img
	 */
	public void negativeTransformation(byte[] img) {
		for (int i=0;i<img.length;i++)
			img[i] = (byte)(255 - (int)(img[i] & 0xFF));
	}
	/**
	 * Task 3:
     * Implement the log transformation.
     * @param img
     */
    public void logTransformation(byte[] img) {
    	double c = 255 / Math.log(256);
		byte[] tmp = new byte[256];
		for (int i=0;i<256;i++)
			tmp[i] = (byte)(c*Math.log(i+1));
		for (int i=0;i<img.length;i++)
			img[i] = tmp[(int)(img[i] & 0xFF)];
    }
    /**
     * Task 4:
     * Implement bit-plane slicing
     *
     * @param img
     * @param mask - between 0 to 255 in decimal
     */
	public void bitPlaneSlicing(byte img[], int mask) {
        System.out.println("Using mask: " + Integer.toBinaryString(mask));
        int tmp[] = new int[256];
        for (int i=0;i<tmp.length;i++)
        	tmp[i] = (int)((i & mask) / (double)mask * 255);
     	for (int i=0;i<img.length;i++)
     		img[i] = (byte)tmp[img[i] & 0xFF];
    }
	/**
	 * Task 5:
	 * Calculate the histogram of the image.
	 * @param img
	 * @return the histogram
	 */
	public int[] histogram (byte[] img) {
		int[] returnValue = new int[256];
		for (int i=0;i<img.length;i++)
			returnValue[img[i]&0xFF]++;
		return returnValue;
	}

	/**
	 * Task 6:
	 * Calculate the cumulative histogram of the image.
	 * @param img
	 * @return the histogram
	 */
	public int[] cumulativeHistogram (byte[] img) {
		int[] histogram = histogram(img);
		int[] returnValue = new int[256];
		for (int i=0;i<histogram.length;i++)
			for (int j=0;j<i;j++) 
				returnValue[i] += histogram[j];
		return returnValue;
	}

	/**
	 * Task 7:
	 * Perform histogram equalization.
	 * @param img
	 */
	public void histogramEqualization(byte img[]) {
		int[] T = cumulativeHistogram(img);
		for (int i=0;i<T.length;i++) 
			T[i] *= 255.0/img.length;
		for (int i=0;i<img.length;i++)
			img[i] = (byte)T[(int)(img[i] & 0xFF)];
	}
	/**
	 * Homework 1:
	 * Implement the box smoothing filter.
	 *
	 * @param img the graylevel image (row major representation)
	 * @param w width of the image
	 * @param h height of the image
	 * @param filterSize the size of the filter, which is supplied by the user
	 */
	public void boxSmoothFilter(byte[] img, int w, int h, int filterSize) {
		System.out.println("TODO: Homework 1");
	}
	/**
	 * Homework 2:
	 * Implement the median filter. 
	 * The java function java.util.Arrays.sort(array) can be used to sort an array.
	 *
	 * @param img the graylevel image (row major representation)
	 * @param w width of the image
	 * @param h height of the image
	 * @param filterSize the size of the filter, which is supplied by the user
	 */
	public void medianFilter(byte[] img, int w, int h, int filterSize) {
		System.out.println("TODO: Homework 2");
	}
	/** 
	 * Homework 3:
	 * Implement the Laplacian filter (isotropic mask for rotations in increments of 45 deg)
	 *
	 * @param img the graylevel image (row major)
	 * @param w width of the image
	 * @param h height of the image
	 */
	public void laplacianFilter(byte img[], int w, int h) {
		System.out.println("TODO: Homework 3");
	}
}
