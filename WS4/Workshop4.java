public class Workshop4 {

   /**
	 * Task 1:
	 * Display the gradient image.
	 *
	 * @param f the graylevel image (row major representation)
	 * @param width of the image
	 * @param height of the image
	 */
	public void gradientImage(byte[] f, int width, int height) {
		double gradient[] = new double[f.length];
		for (int x = 1;x<height-1;x++) {
			for (int y = 1;y<width-1;y++) {
				int fx = ((int)f[(x+1)*width+y]&0xFF) - ((int)f[(x-1)*width+y]& 0xFF);
				int fy = ((int)f[x*width+y+1]&0xFF) - ((int)f[x*width+y-1]& 0xFF);
				gradient[x*width+y] = Math.sqrt(fx*fx + fy*fy);
			}
		}
		double max = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < gradient.length; i++)
			max = Math.max(gradient[i], max);
		double c = 255.0 / max;
		for (int i = 0; i < gradient.length; i++)
			f[i] = (byte)(gradient[i] * c);
	}

    /**
	 * Task 2:
	 * Display the response map R of the harris corner detector. For simplicity, you may use a 5x5 box smoothing filter instead of Gaussian.
	 * Set all pixels to white for which R > 100000.
	 *
	 * @param f the graylevel image (row major representation)
	 * @param width of the image
	 * @param height of the image
	 */
	public void cornerResponseImage(byte[] f, int width, int height) {
		double fx2[] = new double[f.length];
		double fy2[] = new double[f.length];
		double fxy[] = new double[f.length];
		double R[] = new double[f.length];
		for (int x = 1;x<height-1;x++) {
			for (int y = 1;y<width-1;y++) {
				int fx = ((int)f[(x+1)*width+y]&0xFF) - ((int)f[(x-1)*width+y]& 0xFF);
				int fy = ((int)f[x*width+y+1]&0xFF) - ((int)f[x*width+y-1]& 0xFF);
				fx2[x*width+y] = fx * fx;
				fy2[x*width+y] = fy * fy;
				fxy[x*width+y] = fx * fy;
			}
		}

		double[] filter = new double []{1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1};
		for (int i=0;i<filter.length;i++)
			filter[i] /= 25.0;
		convolve(fx2, width, height, filter, 5);
		convolve(fy2, width, height, filter, 5);
		convolve(fxy, width, height, filter, 5);


		for (int x = 1;x<height-1;x++) {
			for (int y = 1;y<width-1;y++)
				R[x*width+y] = fx2[x*width+y] * fy2[x*width+y] - fxy[x*width+y] * fxy[x*width+y] - 0.04
						* ((fx2[x*width+y] + fy2[x*width+y]) * (fx2[x*width+y] + fy2[x*width+y]));	
		}

		for (int i = 0; i < R.length; i++) {
			if (R[i] > 100000) f[i] = (byte)255;
			else f[i] = (byte)0;
		}

	}

	private void convolve(double[] img, int w, int h, double filter[], int filterSize) {
		double[] copy = new double[img.length];
		System.arraycopy(img, 0, copy, 0, img.length);
		int halfFilterSize = filterSize/2;
		for (int y=halfFilterSize;y<h-halfFilterSize;y++){
			for (int x=halfFilterSize;x<w-halfFilterSize;x++) {
				int sum = 0;
				for (int v=-halfFilterSize;v<=+halfFilterSize;v++){
					for (int u=-halfFilterSize;u<=+halfFilterSize;u++)
						sum += filter[(v+halfFilterSize)*filterSize+(u+halfFilterSize)]*(copy[(y+v)*w+(x+u)]);
				}
				img[y*w+x] = sum;
			}
		}
	}


}
