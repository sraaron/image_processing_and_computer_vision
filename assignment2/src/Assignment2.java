import java.util.ArrayList;
/**
 * Assignment 2
 * You may only edit and submit this file to moodle before the deadline.
 * NAME: AARON, Sheshan Ryan
 */
public class Assignment2 {
	/**
	 * Task 1
	 * Implement 2D Gaussian smoothing. Hints:
	 * - Perform 2D filtering by applying 1D filters.
	 * - Compute and use a proper filter size for a 1D Gaussian mask based on the sigma parameter.
	 * - Employ partial filters to handle the boundary of the image
	 * 5 Marks
	 */
	public void gaussianSmooth(final byte[] f, int width, int height, double sigma) {
		double[] img = new double[width * height];
		byteToDouble(f, img);
		gaussianSmooth(img, width, height, sigma);
		doubleToByte(f, img);
		//System.out.println("finish");
	}

	private static void byteToDouble(byte[] f, double[] img){
		for(int i=0; i<img.length; i++)
			img[i] = (double)(1.0 * (f[i] & 0xff));
	}

	private static void doubleToByte(byte[] f, double [] img){
		//rescaling values
		double max = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < img.length; i++)
			max = Math.max(img[i], max);
		double c = 255.0 / max;
		for (int i = 0; i < img.length; i++)
			f[i] = (byte)(c*img[i]);
	}

	private static void gaussianSmooth(double[] img, int width, int height, double sigma) {
		//create gaussian filter
		//calculate radius
		int radius = (int)Math.sqrt((6.908)*(2*sigma*sigma));
		// compute coefficients of the gaussian filter
		double[] filter = new double[2*radius+1];
		double filtersum = 0;

		for(int i=-radius;i<=radius;i++) {
			double g = gaussian(i,sigma);
			filter[i+radius]=g;
			filtersum+=g;
		}

		for(int i=-radius;i<=radius;i++) {
			filter[i+radius]=filter[i+radius]/(2*filtersum);//2 since 1D-filter will be used in x and y directions
		}

		convolveInX(img, width, height, filter, radius);
		convolveInY(img, width, height, filter, radius);

	}

	private static void convolveInX(double[] img, int width, int height, double[] filter, int radius){
		double[] copy = new double[img.length];
		System.arraycopy(img, 0, copy, 0, img.length);
		for (int y=0; y<height; y++) {
			for (int x=0; x<width; x++) {
				double rv = 0.0;
				for(int dx=-radius;dx<=radius;dx++) {
					int xk = x + dx;
					//partial filter handling
					if(xk<0)
						xk = -xk;
					if(xk>=width)
						xk = x - dx;
					rv+=filter[dx+radius]*copy[y*width + xk];
				}
				img[y*width+x]  = rv;
			}
		}
	}

	private static void convolveInY(double[] img, int width, int height, double[] filter, int radius){
		double[] copy = new double[img.length];
		System.arraycopy(img, 0, copy, 0, img.length);
		for (int x=0; x<width; x++) {
			for (int y=0; y<height; y++) {
				double rv = 0.0;
				for(int dy=-radius;dy<=radius;dy++) {
					int yk = y + dy;
					//partial filter handling
					if(yk<0)
						yk = -yk;
					if(yk>=height)
						yk = y - dy;
					rv+=filter[dy+radius]*copy[x+width*yk];
				}
				img[y*width+x]  = rv;
			}
		}
	}

	private static double gaussian(double x, double sigma) {
		double exp = (x*x)/(2*sigma*sigma);
		double coeff = 1.0/(2*Math.PI*sigma*sigma);
		double rv = coeff*Math.exp( -exp);
		return rv;
	}
	/**
	 * Task 2: 
	 * Implement the Harris Corner Detection Algorithm. Hints:
	 * - Make proper use of the sigma and threshold values.
	 * - Suppress the non-maximal corner candidates
	 * - Use Gaussian smoothing on the squared derivative images
	 * - Computer corners up to sub-pixel accuracy
	 * 5 Marks
	 */
	public ArrayList<double[]> detectCorners(final byte[] f, int width, int height, double sigma, double threshold) {
		ArrayList<double[]> cornersOut = new ArrayList<double[]>();
		/**
		 * Parameter values that worked for both images in detecting corners of the black squares:
		 * signma = 1.5
		 * threshold = 700,000
		 */
		//double[] aTestCorner = new double[2];
		//double x = 10.0;
		//double y = 100.0;
		//aTestCorner[0] = x;
		//aTestCorner[1] = y;
		//cornersOut.add(aTestCorner);

		double fx2[] = new double[f.length];
		double fy2[] = new double[f.length];
		double fxy[] = new double[f.length];
		double R[] = new double[f.length];
		for (int y = 1;y<height-1;y++) {
			for (int x = 1;x<width-1;x++) {
				int fy = ((int)f[(y+1)*width+x]&0xFF) - ((int)f[(y-1)*width+x]& 0xFF);
				int fx = ((int)f[y*width+x+1]&0xFF) - ((int)f[y*width+x-1]& 0xFF);
				fx2[y*width+x] = fx * fx;
				fy2[y*width+x] = fy * fy;
				fxy[y*width+x] = fx * fy;
			}
		}

		gaussianSmooth(fx2, width, height, sigma);
		gaussianSmooth(fy2, width, height, sigma);
		gaussianSmooth(fxy, width, height, sigma);

		R = computeR(fx2, fy2, fxy, width, height);
		// supress non-maximal corners
		for (int y=1; y<height-1; y++) {
			for (int x=1; x<width-1; x++) {
				double hh = R[x+width*y];
				if (hh<threshold) continue;
				if (!isLocalMaxima(R, x, y, width)) continue;

				// Sub-pixel Accuracy
				double y2 = y - (R[x+width*(y+1)] - R[x+width*(y-1)]) / (2 * (R[x+width*(y+1)] + R[x+width*(y-1)] - 2 * R[x+width*(y)]));
				double x2 = x - (R[(x+1)+width*y] - R[(x-1)+width*y]) / (2 * (R[(x+1)+width*y] + R[(x-1)+width*y] - 2 * R[(x)+width*y]));

				double[] aTestCorner = new double[2];
				//double x = 10.0;
				//double y = 100.0;
				aTestCorner[0] = x2;
				aTestCorner[1] = y2;
				cornersOut.add(aTestCorner);
			}
		}

		return cornersOut;
	}

	public double[] computeR(double[] fx2, double[] fy2, double[] fxy, int width, int height)
	{
		double k = 0.04;
		double[] rv = new double[width*height];
		for (int y=0; y<height; y++) {
			for (int x=0; x<width; x++) {
				double m00 = fx2[y*width+x];
				double m01 = fxy[y*width+x];
				double m10 = fxy[y*width+x];
				double m11 = fy2[y*width+x];
				rv[y*width+x] = m00*m11 - m01*m10 - k*(m00+m11)*(m00+m11);
			}
		}
		return rv;
	}

	private boolean isLocalMaxima(double[] R, int x, int y, int width) {
		int n=8;
		int[] dx = new int[] {-1,0,1,1,1,0,-1,-1};
		int[] dy = new int[] {-1,-1,-1,0,1,1,1,0};
		double w =  R[x + width*y];
		for(int i=0;i<n;i++) {
			double wk = R[(x+dx[i])+ width*(y+dy[i])];
			if (wk>=w) return false;
		}
		return true;
	}

	/**
	 * This method is called when the Task 3 is clicked in the UI.
	 * This method together with various helper methods will generate many 2D to 3D correspondences.
	 * No need to edit this function.
	 */
	public Matrix determineProjectionMatrix(ArrayList<double[]> corners, ArrayList<double[]> pointList2D, ArrayList<double[]> pointList3D) {
		Matrix projectionMatrix = new Matrix();
		ArrayList<double[]> yzPlane2DPoints = new ArrayList<double[]>();
		ArrayList<double[]> yzPlane3DPoints = new ArrayList<double[]>();
		ArrayList<double[]> xzPlane2DPoints = new ArrayList<double[]>();
		ArrayList<double[]> xzPlane3DPoints = new ArrayList<double[]>();
		for (int i=0;i<pointList3D.size();i++) {
			if (pointList3D.get(i)[1] == 0) {
				xzPlane2DPoints.add(pointList2D.get(i));
				xzPlane3DPoints.add(pointList3D.get(i));
			} else if (pointList3D.get(i)[0] == 0) {
				yzPlane2DPoints.add(pointList2D.get(i));
				yzPlane3DPoints.add(pointList3D.get(i));
			}
		}
		ArrayList<double[]> final2D = new ArrayList<double[]>();
		ArrayList<double[]> final3D = new ArrayList<double[]>();
		if ((xzPlane2DPoints.size() >= 4 && yzPlane2DPoints.size() >= 4)) {
			double gridSize = 0.5;
			double max_u, max_v;
			int[] gridCount = new int[3];
			if (yzPlane2DPoints.size() >= 4) {
				Matrix P_YZ = performPlanarCalibration(yzPlane2DPoints, yzPlane3DPoints, 0);
				max_u = 0;
				max_v = 0;
				for (int i = 0; i < yzPlane3DPoints.size(); i++) {
					max_u = Math.max(yzPlane3DPoints.get(i)[1], max_u);
					max_v = Math.max(yzPlane3DPoints.get(i)[2], max_v);
				}
				gridCount[0] = 1;
				gridCount[1] = (int) Math.floor(max_u / gridSize + 1.5);
				gridCount[2] = (int) Math.floor(max_v / gridSize + 1.5);
				obtainRefiningPoints(corners, P_YZ, gridSize, gridCount, final2D, final3D);
			}
			if (xzPlane2DPoints.size() >= 4) {
				Matrix P_XZ = performPlanarCalibration(xzPlane2DPoints, xzPlane3DPoints, 1);
				max_u = 0;
				max_v = 0;
				for (int i = 0; i < (int) xzPlane3DPoints.size(); i++) {
					max_u = Math.max(xzPlane3DPoints.get(i)[0], max_u);
					max_v = Math.max(xzPlane3DPoints.get(i)[2], max_v);
				}
				gridCount[0] = (int) Math.floor(max_u / gridSize + 1.5);
				gridCount[1] = 1;
				gridCount[2] = (int) Math.floor(max_v / gridSize + 1.5);
				obtainRefiningPoints(corners, P_XZ, gridSize, gridCount, final2D, final3D);
			}
		}
		corners.clear();
		corners.addAll(final2D);
		pointList2D.clear();
		pointList2D.addAll(final2D);
		pointList3D.clear();
		pointList3D.addAll(final3D);
		return performCalibration(final2D, final3D);
	}

	/**
	 * No need to edit this function.
	 */
	private Matrix performPlanarCalibration(ArrayList<double[]> points2d_in, ArrayList<double[]> points3d_in, int planeID) {
		int numberOfPoints = points2d_in.size();
		Matrix A = new Matrix(numberOfPoints * 2, 8);
		Matrix B = new Matrix(numberOfPoints * 2, 1);
		int coordinateIndex = 0;
		if (planeID == 0) coordinateIndex = 1;
		if (planeID == 1) coordinateIndex = 0;
		for (int i = 0; i < numberOfPoints; i++) {
			int c = 0;
			A.set(2*i, c++, points3d_in.get(i)[coordinateIndex]);
			A.set(2*i, c++, points3d_in.get(i)[2]);
			A.set(2*i, c++, 1.0);
			A.set(2*i, c++, 0.0);
			A.set(2*i, c++, 0.0);
			A.set(2*i, c++, 0.0);
			A.set(i << 1, c++, -points2d_in.get(i)[0] * points3d_in.get(i)[coordinateIndex]);
			A.set(i << 1, c++, -points2d_in.get(i)[0] * points3d_in.get(i)[2]);
			c = 0;
			A.set(i * 2 + 1, c++, 0.0);
			A.set(i * 2 + 1, c++, 0.0);
			A.set(i * 2 + 1, c++, 0.0);
			A.set(i * 2 + 1, c++, points3d_in.get(i)[coordinateIndex]);
			A.set(i * 2 + 1, c++, points3d_in.get(i)[2]);
			A.set(i * 2 + 1, c++, 1.0);
			A.set(i * 2 + 1, c++, -(points2d_in).get(i)[1] * (points3d_in).get(i)[coordinateIndex]);
			A.set(i * 2 + 1, c++, -(points2d_in).get(i)[1] * (points3d_in).get(i)[2]);
			B.set(i * 2, points2d_in.get(i)[0]);
			B.set(i * 2 + 1, points2d_in.get(i)[1]);
		}
		Matrix x = A.inverse().mul(B);
		Matrix prj = new Matrix(3, 3);
		for (int i = 0; i < 8; i++) {
			prj.set(i / 3, i % 3, x.get(i, 0));
		}
		prj.set(2, 2, 1.0);
		Matrix P = new Matrix(3, 4);
		for (int i = 0; i < 3; i++) {
			int c = 0;
			for (int j = 0; j < 4; j++) {
				if (j == planeID)
					P.set(i, j, 0.0);
				else
					P.set(i, j, prj.get(i, c++));
			}
		}
		return P;
	}
	/**
	 * No need to edit this function.
	 */
	private void obtainRefiningPoints(ArrayList<double[]> points2d_in, Matrix P, double gridSize_in, int gridCount_in[], ArrayList<double[]> points2d_out, ArrayList<double[]> points3d_out) {
		double fuzziness = 25;
		for (int x = (gridCount_in[0] > 1 ? 1 : 0); x < gridCount_in[0]; x += 2) {
			for (int y = (gridCount_in[1] > 1 ? 1 : 0); y < gridCount_in[1]; y += 2) {
				for (int z = 1; z < gridCount_in[2]; z += 2) {
					Matrix x3D = new Matrix(4, 1);
					x3D.set(0, 0, x * gridSize_in);
					x3D.set(1, 0, y * gridSize_in);
					x3D.set(2, 0, z * gridSize_in);
					x3D.set(3, 0, 1.0);
					Matrix x2D = P.mul(x3D);
					x2D = x2D.div(x2D.get(2));
					double minDistance = 100000;
					int minDistanceIndex = 0;
					for (int j=0;j<points2d_in.size();j++) {
						double distance = Math.pow(x2D.get(0) - points2d_in.get(j)[0], 2) + Math.pow(x2D.get(1) - points2d_in.get(j)[1], 2);
						if (distance < minDistance) {
							minDistance = distance;
							minDistanceIndex = j;
						}
					}
					if (minDistance < fuzziness) {
						points2d_out.add(points2d_in.get(minDistanceIndex));
						points3d_out.add(new Matrix(x3D).data);
					}
				}
			}
		}
	}

	/**
	 * Task 3:
	 * Perform camera calibration and return the 3x4 camera projection matrix from the provided
	 * correspondences. Hints:
	 * - Solve the equation Ap = 0 on slide 5 of chapter 7.
	 * - The matrix A is of size (points2d.size() * 2, 12)
	 * - Use the provided function SVD2 in the class Matrix.
	 * - The solution to Ap = 0 will be the last column in V return by SVD2.
	 * 5 Marks
	 */
	private Matrix performCalibration(ArrayList<double[]> points2d, ArrayList<double[]> points3d) {

		System.out.println("out");
		//create A, Chapter 7, pg 5
		double [] m = new double[2*points3d.size()*12];

		int width = 12;
		for(int y = 0; y<points3d.size()*2; y++){
			for(int x=0; x<width; x++){
				if((y%2)==0){
					if(x>=0 && x<4) {
						/*System.out.println("X " + x);
						System.out.println("Y " + y);
						System.out.println("SZ " + points3d.size());
						if(x == 0 && y ==24){
							System.out.println("BREAK");
						}*/
						m[y * width + x] = points3d.get(y/2)[x];
					}
					else if(x>3 && x<8)
						m[y*width+x] = 0;
					else if(x>7 && x<12)
						m[y*width+x] = -(points3d.get(y/2)[x-8] * points2d.get(y/2)[0]);
				}
				else{
					if(x>=0 && x<4)
						m[y*width+x] = 0;
					else if(x>3 && x<8)
						m[y*width+x] = points3d.get(y/2)[x-4];
					else if(x>7 && x<12)
						m[y*width+x] = -(points3d.get(y/2)[x-8] * points2d.get(y/2)[1]);
				}
			}
		}
		//System.out.println("Finish");
		//create matrix A
		Matrix A = new Matrix(points2d.size() * 2, 12, m);//rows,columns

		Matrix U = new Matrix ();
		Matrix D = new Matrix ();
		Matrix V = new Matrix ();

		//perform SVD2 on A
		A.SVD2(U,D,V);

		double[] rv = new double[12];
		Matrix p = new Matrix();
		V.getCol(V.cols-1,p);

		for(int i=0; i<12; i++){
			rv[i] = p.get(i,0)/p.get(11,0);
		}

		return new Matrix(3, 4, rv);
	}
}
