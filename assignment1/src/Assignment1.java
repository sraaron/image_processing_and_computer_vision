/**
  * Assignment 1
  * You may only edit and submit this file to moodle before the deadline.
  * NAME: AARON, Sheshan Ryan
  */
public class Assignment1 {

    private void logTransformation(byte[] img) {
      double c = 255 / Math.log(256);
      byte[] tmp = new byte[256];
      for (int i=0;i<256;i++)
        tmp[i] = (byte)(c*Math.log(i+1));
      for (int i=0;i<img.length;i++)
        img[i] = tmp[(int)(img[i] & 0xFF)];
    }

  /**
    * Task 1
    * Implement the Fast Fourier Transform (FFT) and display the Fourier spectrum.
    * Results should be equivalent to those obtained in the Workshop 3.
    * The implementation details of the FFT can be obtained from section 4.11 of chapter 4.
    * 6 Marks
    */
    public void fourierTransform(byte[] img, int width, int height) {
      Complex[] f = imgToComplex(img, width, height);
      center(f, width, height);
      Complex[] F = fft2D(f, width, height);
      double max = Double.NEGATIVE_INFINITY;
      for (int i = 0; i < F.length; i++)
        max = Math.max(F[i].getNorm(), max);
      double c = 255.0 / max;
      for (int i = 0; i < img.length; i++)
        img[i] = (byte)(c*F[i].getNorm());

      logTransformation(img);
      /*for(int x =0; x<height; x++){
        for(int y = 0; y<width; y++){
          System.out.print(x + "," + y + ",");
          System.out.println(img[x * width + y] & 0xFF);
        }
      }*/
    }

    private Complex[] imgToComplex(byte[] img, int width, int height){
    Complex[] F = new Complex[width*height];
    for(int i =0; i<width*height; i++)
      F[i] = Complex.fromPolar((double) (img[i] & 0xFF), 0);
    return  F;
  }

    private void center(Complex[] f, int width, int height) {
      for (int x = 0; x < height; x++)
        for (int y = 0; y < width; y++)
          f[x*width+y].setReal(f[x*width+y].getReal() * Math.pow(-1, (x + y)));
    }

    private Complex[] ft1D(Complex[] fx, int width) {
      Complex[] F = new Complex[width];
      for (int v = 0; v < width; v++) {
        F[v] = new Complex(); //create new Complex object
          for (int y = 0; y < width; y++) {
            //double a = (double) (img[y] & 0xFF) * Math.pow(-1, y);
            double b = -2.0 * Math.PI * (v * y / (double) width);
            Complex c = Complex.fromPolar(1, b).mul(fx[y]);
            F[v] = F[v].plus(c);
          }
        }
      return F;
    }

    private Complex[] sepft(Complex[] cmplxImg, int width, int height) {
      Complex[] F = new Complex[width*height];
      Complex[] fx = new Complex[width];
      Complex[] Fv;
      for (int x = 0; x < height; x++) {
        for (int y = 0; y < width; y++)
          fx[y] = cmplxImg[x * width + y]; //each row
        Fv = ft1D(fx, width);
        for (int y = 0; y < width; y++)
          F[x * width + y] = Fv[y];
      }
      Complex[] fy = new Complex[height];
      Complex[] Fu;
      for (int y = 0; y < width; y++) {
        for (int x = 0; x < height; x++)
          fy[x] = F[x * width + y]; //each column
        Fu = ft1D(fy, width);
        for (int x = 0; x < height; x++)
          F[x * width + y] = Fu[x];
      }
      return F;
    }

    private Complex[] ft(byte[] img, int width, int height) {
      Complex[] F = new Complex[width*height];
      for (int u = 0;u<height;u++) {
        for (int v = 0;v<width;v++) {
          F[u*width + v] = new Complex();
          for (int x = 0;x<height;x++) {
            for (int y = 0;y<width;y++) {
              double a = (double)(img[x*width+y]&0xFF)*Math.pow(-1, x+y);
              double b = -2.0 * Math.PI * (u*x/(double)height + v*y/(double)width);
              Complex c = Complex.fromPolar(1, b).mul(a);
              F[u*width+v] = F[u*width+v].plus(c);
            }
          }
        }
      }
    return F;
  }

    private void fft1D(Complex[] f1D, Complex[] rvF, int twoK) {
      // base
      if (twoK <= 1) {
        rvF[0] = f1D[0];
        return;
      }
      Complex[] oddF = new Complex[twoK / 2];
      for (int k = 0; k < twoK / 2; k++)
        oddF[k] = f1D[2 * k + 1];

      Complex[] evenF = new Complex[twoK / 2];
      for (int k = 0; k < twoK / 2; k++)
        evenF[k] = f1D[2 * k];

      Complex[] rvEven = new Complex[twoK / 2];
      fft1D(evenF, rvEven, twoK / 2);

      Complex[] rvOdd = new Complex[twoK / 2];
      fft1D(oddF, rvOdd, twoK / 2);

      for (int k = 0; k < twoK / 2; k++) {
        double kpi = -2 * k * Math.PI / twoK;
        Complex wk = new Complex(Math.cos(kpi), Math.sin(kpi));
        rvF[k] = rvEven[k].plus(wk.mul(rvOdd[k]));
        rvF[k + twoK / 2] = rvEven[k].minus(wk.mul(rvOdd[k]));
      }

    }

    private Complex[] fft2D(Complex[] cmplxImg, int width, int height) {
      Complex[] F = new Complex[width*height];
      Complex[] fx = new Complex[width];
      Complex[] Fv = new Complex[width];
      for (int x = 0; x < height; x++) {
        for (int y = 0; y < width; y++)
          fx[y] = cmplxImg[x * width + y]; //each row

        fft1D(fx, Fv, width);
        for (int y = 0; y < width; y++)
          F[x * width + y] = Fv[y];
          //System.out.println(Fv[y]);

      }
      Complex[] fy = new Complex[height];
      Complex[] Fu = new Complex[height];
      for (int y = 0; y < width; y++) {
        for (int x = 0; x < height; x++)
          fy[x] = F[x * width + y]; //each column

        fft1D(fy, Fu, height);
        for (int x = 0; x < height; x++)
          F[x * width + y] = Fu[x];
          //System.out.println(Fu[x]);
      }
      return F;
    }    

    /**
      * Task 2
      * Describe the appearance of the Fourier spectrum that is obtain from the image created in this function.
      * Why does the Fourier spectrum appear as it does? Provide a detailed explanation as a comment to this function.
      * No need to implement any code here.
      * 2 Marks
      */
  public void changeImage(byte[] img, int width, int height) {
    for (int x = 0;x<height;x++)
      for (int y = 0;y<width;y++)
        img[x*width+y] = (byte)(255.0 * Math.cos(x / 4.0 * Math.PI));
    /*
    The spectrum obtained from the image created in this function consists of 8 evenly spaced dots
    going vertically through the center of the image. The image size is 128 x 128 pixels.

    The image that is created is a cosine function of pi/4 which travels down along the x axis.
    The brightest region represents the peaks, while the darkest region represents the troughs.

    The image function can be represented in euler form as:
      0.5e^-(0.25*j*pi*x) + 0.5e^(0.25*j*pi*x)

    The Fourier Transform decomposes the image into sinusoids.
    The 2D Fourier Transform can be thought of as getting the Fourier Transform of each row in the
    image and then getting the Fourier Transform of each column of the result.

    Each row in the image is a constant value, hence the Fourier Transform of each row of the image
    would be just a single frequency at the center.

    Therefore the result of the 1D FT would just be values of cos(x / 4.0 * PI) on the central column.

    DTFT is perdiodic with period 2*PI. The whole image plane represents 2*PI x 2*PI angular distance.
    The sampling frequency (omega_s) of the Fourier Transform is:
      omega_s = 2*PI / Ts, where Ts, the period = 2*PI
    so omega_s = 1

    The frequency of the cosine function in the image would then be:
      omega_img = 2 * PI / (PI / 4) = 8
    The image has higher frequency than the Fourier Transform and hence by the Nyquist-Shannon Sampling
    theorem we can see that this would cause aliasing. Hence why we see 8 dots vertically along the central
    column of the image's Fourier Transform. Since here, the Fourier Transform is sampled at at PI/4, so we
    get an impulse (since FT of cosine is an impulse centered 0) and shifted PI/4 multiples over 2*PI. Which
    has 8 PI/4 in 2*PI.
     */

    /*
    Initial ideas: Sampling frequency vs peridocity of lines (some sort of aliasing)
    FT will get the change, change is along vertical axis
     */
  }

    /**
      * Task 3
      * Perform second order (n=2) ButterWorth low pass filtering in the frequency domain.
      * 5 Marks
      */
  public void filtering(byte[] img, int width, int height, double d0) {
    Complex[] F;
    Complex [] cmplxImg = imgToComplex(img, width, height);
    center(cmplxImg, width, height); // step 1
    F = fft2D(cmplxImg, width, height); // step 2
    blpf(F, width, height, d0); // step 3
    ifft2D(F, cmplxImg, width, height);// step 4 and 5
    center(cmplxImg, width, height); // step 6
    for(int i =0; i<width*height; i++)
        img[i] = (byte)cmplxImg[i].getReal();
  }

  //BLPF second order (n=2)
  private void blpf(Complex[] F, int width, int height, double D0) {
    int cY = width / 2;
    int cX = height / 2;
    //int filterOrder = 2;
    for (int y = 0; y < width; y++)
      for (int x = 0; x < height; x++){
        double temp;
        //double euclideanDist = Math.pow(Math.pow(x-cX, 2) + Math.pow(y-cY, 2), filterOrder/2);
        //temp = euclideanDist/Math.pow(D0, filterOrder);
        double euclideanDist = Math.pow(x-cX, 2) + Math.pow(y-cY, 2);
        temp = euclideanDist/Math.pow(D0, 2);
        temp = 1 + temp;
        temp = 1/temp;
        F[x*width+y] = F[x*width+y].mul(temp);
      }
  }

  private void ifft2D(Complex[] F, Complex[] f, int width, int height) {

    Complex[] FStar = new Complex[width*height];
    Complex[] rv;

    for (int i = 0; i < width*height; i++)
      FStar[i] = F[i].getConjugate();

    rv = fft2D(FStar, width, height);
    for (int i = 0; i < width*height; i++){
      f[i] = rv[i].getConjugate();
      f[i] = f[i].div((double)width*height);// 1/MN
    }
  }

  //Butterworth notch reject filter
  private void bhpf(Complex[] F, int width, int height, double D0, int uk, int vk) {
    int cY = width / 2;
    int cX = height / 2;
    //int filterOrder = 2;
    for (int y = 0; y < width; y++)
      for (int x = 0; x < height; x++){
        double val1 = 0;
        double val2 = 0;
        double value = 0;
        //double Dk = Math.pow(Math.pow(x-cX-uk, 2) + Math.pow(y-cY-vk, 2),filterOrder/2);
        //double D_k = Math.pow(Math.pow(x-cX+uk, 2) + Math.pow(y-cY+vk, 2),filterOrder/2);
        //val1 = Math.pow(D0,filterOrder)/Dk;
        //val2 = Math.pow(D0,filterOrder)/D_k;
        double Dk = Math.pow(x-cX-uk, 2) + Math.pow(y-cY-vk, 2);
        double D_k = Math.pow(x-cX+uk, 2) + Math.pow(y-cY+vk, 2);
        val1 = Math.pow(D0,2)/Dk;
        val2 = Math.pow(D0,2)/D_k;
        val1 = 1 + val1;
        val2 = 1 + val2;
        val1 = 1/val1;
        val2 = 1/val2;
        value = val1 * val2;
        F[x*width+y] = F[x*width+y].mul(value);
      }
  }

  /*
  //returns median. Assumes odd number of elements in array
  private byte getMedian(byte[] byteArray){
    Arrays.sort(byteArray);
    return byteArray[byteArray.length/2 +1];
  }

  private byte[] medianFilter(byte[] img, int width, int height) {
    byte[] rvImg = new byte[width*height];
    int N = 2;
    byte[] window = new byte[(2*N + 1) * (2*N + 1)];
    for(int y=N; y<width-N; y++){
      for(int x=N; x<height-N; x++){
        int i =0;
        for(int ny = -N; ny<N+1; ny++) {
          for (int nx = -N; nx <N+1; nx++) {
            window[i] = img[(x+nx)*width + (y+ny)];
            i++;
          }
        }
        rvImg[x*width + y] = getMedian(window);
      }
    }
    return  rvImg;
  }
  */
    /**
      * Task 4
      * Apply a suitable filter to the car image to attenuate the "impulse-like" bursts in the image.
      * 2 Mark
      */
  public void filtering2(byte[] img, int width, int height) {
    //Section 4.10.2 in DIP textbook discusses a technique for reducing moire interference patterns
    //Filter using a Notch Filter
    Complex[] F;
    Complex [] cmplxImg = imgToComplex(img, width, height);
    center(cmplxImg, width, height); // step 1
    F = fft2D(cmplxImg, width, height); // step 2
    //do a Low Pass Filter first to remove any high frequency noise
    blpf(F, width, height, 40); // step 3
    //for(int q = 0; q<1; q++){
      bhpf(F, width, height, 10, 20, 22); // step 3
      bhpf(F, width, height, 10, -22, 21); // step 3
      bhpf(F, width, height, 10, 0, 41); // step 3
    //bhpf(F, width, height, 3, 20, 85); // step 3
    //bhpf(F, width, height, 3, 23, 41); // step 3
    //}
    //f = F;
    ifft2D(F, cmplxImg, width, height);// step 4 and 5
    center(cmplxImg, width, height); // step 6

    double max = Double.NEGATIVE_INFINITY;
    for (int i = 0; i < cmplxImg.length; i++)
      max = Math.max(cmplxImg[i].getReal(), max);
    double c = 255.0 / max;
    for (int i = 0; i < img.length; i++)
      img[i] = (byte)(c*cmplxImg[i].getReal());
    //byte[] rvImg = img;
    //int N = 0;
    //img = medianFilter(img, width, height);
    /*for(int y=N; y<width-N; y++) {
      for (int x = N; x < height - N; x++) {
        System.out.println(rvImg[x*width + y] - img[x*width + y]);
      }
    }*/
    //logTransformation(img);
  }

}
