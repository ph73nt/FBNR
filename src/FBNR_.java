package FBNR;

///////////////////////////////////////////////////////////////////////////////
// FBNR - Fourier Block Noise Reduction
///////////////////////////////////////////////////////////////////////////////
// 
// Neil Thomson, March 2010
// Dept Medical Physics, Nuc Med section,
// Kent and Canterbury Hospital, CT1 3NG, UK.
//
// This plugin implements an FBNR filter on the top image.  The plugin is based
// on Matthew J Guy's article in Nuclear Medicine Communications:
//
//   Guy, M. J., Nuc Med Comm, Fourier block noise reduction: an adaptive filter
//    for reducing Poisson noise in scintigraphic images, March 2008, vol 29,
//    issue 3, pp291-297.
//
///////////////////////////////////////////////////////////////////////////////
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.filter.*;
import ij.plugin.frame.*;
import ij.plugin.*;
import java.io.*;
import javax.swing.Timer;
import java.awt.event.*;

public class FBNR_ implements PlugInFilter {

  ////////////////////////////////////////////////////////////////
  // Global declarations
  ////////////////////////////////////////////////////////////////
  // Images
  private ImagePlus imp, imBlock, imResult;
  // Processors
  private ImageProcessor ipBlock, ipResult;
  // Logging
  private boolean logging, canUpdate;
  // Look up tables for referencing rows and columns in the fourier transform
  private int[] refRow, refCol;
  // Set up a rectangle for the fourier transform
  private Rectangle rect = new Rectangle();

  public int setup(String arg, ImagePlus imp) {
    // Convert image to 32 bit float type before assigning processor
    new ImageConverter(imp).convertToGray32();

    this.imp = imp;
    return DOES_8G + DOES_16 + DOES_32;
  }

  public void run(ImageProcessor ip) {
    // Indices
    int i, startI, startJ, j, m, n, p;
    // The length of side of the processing block
    int blockSide;
    // Maximum number of iterations before meltdown
    int maxIterations;
    // Floats to hold pixel values
    float[] pixFilter, pixFilterPrev, pixBlock, pixResult;
    // Catch errors for reporting
    boolean FBNR_error = false, FBNR_warning = false;

    // get pixel array of the main image
    float[] pixIm = (float[]) ip.getPixels();

    // Set up the variance class (like C struct)
    variance var = new variance();

    // Collect some options from the dialogue box (or args)
    GenericDialog FBNR_Opts = doDialogue();
    FBNR_Opts.showDialog();
    if (!FBNR_Opts.wasCanceled()) {

      // Display a progress bar.  The imagej progress bar will not function
      //  correctly when there are multiple calls to it.  In this plugin
      //  the IJ progress bar is called for every FFT, so calling it to
      //  update the progress of this plugin would result in incorrect display
      timer.start();  // Only update every couple of seconds to save CPU time
      String progressFile = progressReport();

      // Set block size
      if (FBNR_Opts.getNextChoice() == "4x4") {
        blockSide = 4;
      } else {
        blockSide = 8;
      }
      int blockSize = (int) Math.pow(blockSide, 2);

      // Set max iterations per block
      maxIterations = (int) FBNR_Opts.getNextNumber();

      // How much will the filter change by when the rate is discovered
      // to be too high?
      float changeRate = (float) FBNR_Opts.getNextNumber();

      // Choose whether to show log output (very slow)
      logging = FBNR_Opts.getNextBoolean();

      //Set up the reference tables
      refRow = new int[(int) Math.pow(2 * blockSide, 2)];
      refCol = new int[(int) Math.pow(2 * blockSide, 2)];
      setUpLookUpTables(blockSide);

      // Setup the block, result and the filter images
      pixBlock = makeBlockImage(blockSide);
      pixResult = makeResultImage();
      pixFilter = makeFilterArray(blockSide);
      pixFilterPrev = new float[4 * blockSize];

      i = 125;//12,60,80,48
      j = 61;//50,50,40,60

      // Scan the image a total of blocksize times and take an average
      // after each scan, the "startscan" place needs to be incremented
      for (n = 0; n < blockSide; n++) {
        if (logging) {
          IJ.log(Integer.toString(n));
        }
        startI = -1 * (blockSide - 1 - n);
        for (p = 0; p < blockSide; p++) {
          startJ = -1 * (blockSide - 1 - p);
          // Scan y axis of image
          for (j = startJ; j < imp.getHeight() + startJ; j += blockSide) {
            //Scan x-axis of image
            for (i = startI; i < imp.getWidth() + startI; i += blockSide) {
              if (logging) {
                String message = "n = " + n + ", startI = " +
                        startI + ", i = " + i;
                IJ.log(message);
              }

              // Fill the block with values from the main image
              var = setBlockValues(i, j, blockSide, pixIm, pixBlock, var, true);

              // Estimate the variance of the block for the first time
              var = getBlockVariance(var, true);

              if (logging) {
                IJ.log("---------------------------------------------");
                IJ.log("Block variance = " + Double.toString(var.tot0));
                IJ.log("Block noise = " + Double.toString(var.noise));
                IJ.log("Block max = " + Double.toString(var.max));
                IJ.log("---------------------------------------------");
              }

              float change = (float) 1 / blockSide;
              // Iterate the filter until it is correct to the nearest...
              boolean continu = true;
              if (var.max > 0) {
                if (var.tot0 < var.noise) {
                  // This will never converge! Keep the original values
                  // and display a warning
                  FBNR_warning = true;
                } else {
                  m = 0;       // to track the number of iterations

                  // Prevent the loop going into meltdown by limiting
                  // the number of iterations.
                  while (continu) {
                    // Fill the block with values from the main image
                    setBlockValues(i, j, blockSide, pixIm,
                            pixBlock, var, false);
                    // Save the previous filter
                    System.arraycopy(pixFilter, 0, pixFilterPrev, 0,
                            pixFilter.length);
                    // Edit the filter
                    pixFilter = changeFilterArray(pixFilter, blockSide, change);
                    // Filter the block
                    filterImage(pixFilter, blockSide);
                    //for (l=0; l<16; l++) {
                    //    IJ.log("Filter: " + Float.toString(pixFilter[l])
                    //    + "Previous: " + Float.toString(pixFilterPrev[l]));
                    //}
                    // Recalculate variance
                    var = getBlockVariance(var, false);
                    m++;
                    // Show log
                    if (logging) {
                      IJ.log("Iteration: " + Integer.toString(m) +
                              "  Change = " + Float.toString(change) +
                              " i=" + Integer.toString(i));
                      IJ.log("Block variance = " + Double.toString(var.tot));
                      IJ.log("Block noise = " + Double.toString(var.noise));
                      IJ.log("Block residual = " + Double.toString(var.res));
                    }
                    // When the residual is larger than the noise we have
                    // reached a point where we need to condsider what's to
                    // be done.
                    // Should the loop terminate?
                    if (var.res > var.noise) {
                      // change the change!
                      change = change / changeRate;
                      // Reset the filter to the previous value
                      // (ie where noise > res)
                      System.arraycopy(pixFilterPrev, 0, pixFilter, 0,
                              pixFilterPrev.length);
                    }
                    if (Math.abs(var.res - var.noise) < 0.1) {
                      continu = false;
                    }
                    if (m >= maxIterations) {
                      continu = false;
                      FBNR_error = true;
                      if (logging) {
                        String strError = "No convergance at i=" +
                                i + ", j=" + j;
                        IJ.log(strError);
                      }
                    }
                  }
                }

                // Shove the block back into the image now that the
                // noise has gone!
                pushBlock(i, j, blockSide, pixResult, pixBlock, var);
                // Reset the filter or we'll be in all sorts of bother
                pixFilter = resetFilterArray(pixFilter);
              }
            } // i
          }  // j
          updateProgress(progressFile, (int) 100 * (p + n * blockSide) *
                  imp.getHeight() / (blockSize * imp.getHeight()));
        } //p
      }  //n
      // Finalise the image
      finaliseImage(blockSize);
      if (FBNR_error & !logging) {
        IJ.showMessage("FBNR Error!", "Errors have occurred.\n" +
                "Try logging mode or more iterations");
      }
      if (FBNR_warning) {
        IJ.showMessage("Warning!", "Some areas have not been " +
                "filtered due to input homogeniety");
      }
      canUpdate = true;
      updateProgress(progressFile, 100);

    }
  }

  GenericDialog doDialogue() {
    GenericDialog FBNR_Opts = new GenericDialog("FBNR options");

    String[] blockOpts = new String[2];
    blockOpts[0] = "4x4";
    blockOpts[1] = "8x8";

    FBNR_Opts.addChoice("Block size", blockOpts, "4x4");
    FBNR_Opts.addNumericField("Max iterations per block", 50, 0);
    FBNR_Opts.addNumericField("Rate of change of filter", 5, 0);
    FBNR_Opts.addCheckbox("Enable logging (slow)", false);

    return FBNR_Opts;
  }

  void finaliseImage(int blockSize) {
    // Rescale pixel values to make image an average, not a sum
    ipResult.multiply((float) 1 / blockSize);
    // Remove any pixel values below zero
    ipResult.min(0);

    // Show the image
    imResult.show();

    // Get max value to set contast
    ImageStatistics stat = imResult.getStatistics();
    ipResult.setMinAndMax(0, stat.max);

    imResult.updateAndRepaintWindow();
  }

  // Copy values from the main image into the processing block
  variance setBlockValues(int i, int j, int blockSide, float pixIm[],
          float pixBlock[], variance var, boolean first) {

    // Get pixel values within the block
    for (int l = 0; l < blockSide; l++) {      // y values
      for (int k = 0; k < blockSide; k++) {    // x-values
        // Assign the pixel values to the block image - only values within
        // the image are allowed to be copied, otherwise force to zero
        if (i >= 0 & j >= 0 & i < imp.getWidth() & j < imp.getHeight()) {
          pixBlock[k + l * blockSide] = pixIm[i + k + (j + l) * imp.getWidth()];
        } else {
          pixBlock[k + l * blockSide] = 0;
        }
      } // k
    }   // l

    if (first) {
      // Ramp up the grey levels of low count areas or accentuate differences
      // of homogeneous areas where the noise is greater than the variance
      ImageStatistics stat = imBlock.getStatistics();
      var.max = stat.max;

      // scale image by 100 for low count areas to reduce amount
      // of negative pixels
      if (var.max < 25) {
        var.scale = (double) 100;
        ipBlock.multiply(var.scale);
      } else {
        var.scale = (double) 1;
      }

      //Reload stats - in case it's a low count area - before testing
      stat = imBlock.getStatistics();
      // Scale the image up to 1000 times in an attempt
      // to accentuate differences
      while ((Math.pow(stat.stdDev, 2) < stat.mean) & var.scale < 1000) {
        var.scale = var.scale * 10;
        ipBlock.multiply(10);
        stat = imBlock.getStatistics();
      }

    } else {
      ipBlock.multiply(var.scale);
    }

    return var;
  }

  // Copy values from the processing block into the main image
  void pushBlock(int i, int j, int blockSide, float[] pixResult,
          float[] pixBlock, variance var) {

    // Rescale the image
    ipBlock.multiply((double) (1 / var.scale));

    // Get pixel values within the block (assuming they are
    // within the image proper)
    for (int l = 0; l < blockSide; l++) {      // y values
      for (int k = 0; k < blockSide; k++) {    // x-values
        if (i >= 0 & j >= 0 & i < imResult.getWidth() &
                j < imResult.getHeight()) {
          // Assign the pixel values to the image
          pixResult[i + k + (j + l) * imResult.getWidth()] +=
                  pixBlock[k + l * blockSide];
        }  // else do nothing
      } // k
    }   // l
  }

  // Calculate the variance of a block
  variance getBlockVariance(variance var, boolean first) {

    // Let ImageJ calculate the statistics
    ImageStatistics stat = imBlock.getStatistics();
    // Integer value "first" instructs the function that values should
    // be stored differently
    var.tot = Math.pow(stat.stdDev, 2);
    if (first) {
      // Noise in the sub block is estimated as the mean value
      var.noise = stat.mean;
      var.tot0 = var.tot;
    } else {
      // Residual value is difference between the initial blcok variance
      // and the current. This is meaningless unless tot0 has been explicitly
      // set during the program.
      var.res = var.tot0 - var.tot;
    }

    return var;
  }

  // Makes an image block and returns an array to the pixel values
  float[] makeBlockImage(int blockSide) {

    // Create a new image with stack size of 1, filled with black
    imBlock = NewImage.createFloatImage("Block", blockSide, blockSide, 1, 1);
    // Collect the imageprocessor
    ipBlock = imBlock.getProcessor();
    ipBlock.setProgressBar(null);
    //imBlock.show();
    // Get the pixels of the block
    float[] pixBlock = (float[]) ipBlock.getPixels();
    // Force first block to contain known pixel values
    for (int i = 0; i < 16; i++) {
      pixBlock[i] = i;
    }
    return pixBlock;
  }

  // Makes an image to assign final values to without
  // corrupting the input image
  float[] makeResultImage() {

    // Create a new image with stack size of 1, filled with black
    imResult = NewImage.createFloatImage("Result", imp.getWidth(),
            imp.getHeight(), 1, 1);
    // Collect the imageprocessor
    ipResult = imResult.getProcessor();
    ipResult.setProgressBar(null);
    // Get the pixels of the block
    float[] pixResult = (float[]) ipResult.getPixels();
    return pixResult;
  }

  void filterImage(float[] pixFilter, int blockSide) {
    FHT fht = newFHT(ipBlock);
    ((FHT) fht).transform();
    customFilter(fht, pixFilter, blockSide);
    doInverseTransform(fht, ipBlock);
  }

  float[] makeFilterArray(int blockSide) {
    int tot = (int) Math.pow(2 * blockSide, 2);
    float[] pixFilter = new float[tot];
    for (int i = 0; i < tot; i++) {
      pixFilter[i] = (float) 1;
    }
    return pixFilter;
  }

  float[] resetFilterArray(float[] pixFilter) {
    for (int i = 0; i < pixFilter.length; i++) {
      pixFilter[i] = 1;
    }
    return pixFilter;
  }

  float[] changeFilterArray(float[] pixFilter, int blockSide, float change) {
    // Change the filter by a set amount.  The filter is always a low-pass
    // filter, so remove high frequencies.  The high frequencies are at the
    // edges of the array.
    int i = 0, j = 0, k = 0;
    int FFT_Side = 2 * blockSide;
    int perimSide = FFT_Side - 1;
    float value;

    // The number of loops to make is the same as the blockSide as the
    // FHT array contains blockSize*4 elements and we want to loop for each
    // "square circle" of elements in 2D
    for (j = 0; j < blockSide; j++) {
      // No point looping if there's no change to commit
      if (change > 0) {
        // Calculate the value to assign for each pixel
        // If the pixel value would go below zero, cap
        // it and carry the surplus to the next level
        if (pixFilter[refCol[k] + refRow[k] * FFT_Side] - change < 0) {
          change = Math.abs(pixFilter[refCol[k] +
                  refRow[k] * FFT_Side] - change);
          value = 0;
        } else {
          value = pixFilter[refCol[k] + refRow[k] * FFT_Side] - change;
          change = 0;
        }

        // Set pixel values for the perimiter
        int perimeter = 4 * perimSide;
        for (i = k; i < perimeter + k; i++) {
          pixFilter[refCol[i] + refRow[i] * FFT_Side] = value;
        }
        k = i;
        perimSide = perimSide - 2;
      }
    }

    return pixFilter;
  }

  void doInverseTransform(FHT fht, ImageProcessor ip) {
    fht.inverseTransform();
    fht.resetMinAndMax();
    ImageProcessor ip2 = fht;
    fht.setRoi(rect.x, rect.y, rect.width, rect.height);
    ip2 = fht.crop();
    int bitDepth = fht.originalBitDepth >
            0 ? fht.originalBitDepth : imp.getBitDepth();
    switch (bitDepth) {
      case 8:
        ip2 = ip2.convertToByte(true);
        break;
      case 16:
        ip2 = ip2.convertToShort(true);
        break;
      case 24:
        fht.rgb.setBrightness((FloatProcessor) ip2);
        ip2 = fht.rgb;
        fht.rgb = null;
        break;
      case 32:
        break;
    }
    ip.insert(ip2, 0, 0);
  }

  FHT newFHT(ImageProcessor ip) {
    FHT fht;
    int width = ip.getWidth();
    int height = ip.getHeight();
    int maxN = Math.max(width, height);
    int size = 2;
    while (size < 1.5 * maxN) {
      size *= 2;
    }
    rect.x = (int) Math.round((size - width) / 2.0);
    rect.y = (int) Math.round((size - height) / 2.0);
    rect.width = width;
    rect.height = height;
    FFTFilter fftFilter = new FFTFilter();
    fht = new FHT(fftFilter.tileMirror(ip, size, size, rect.x, rect.y));
    //fht.originalWidth = originalWidth;
    //fht.originalHeight = originalHeight;
    //fht.originalBitDepth = imp.getBitDepth();
    return fht;
  }

  void setUpLookUpTables(int blockSide) {
    // The fast hartley transform returns no imaginary values.
    // The number of values in the FHT array is equal to (2*blockSide)^2.
    // Edge values need to be edited so need to create a look up table
    //   to refer to ponts within the image array.

    if (blockSide == 4) {
      int[] row = {0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 3, 4, 5, 6, 7,
        7, 7, 7, 7, 7, 7, 7,
        6, 5, 4, 3, 2, 1,
        1, 1, 1, 1, 1, 1,
        2, 3, 4, 5, 6,
        6, 6, 6, 6, 6,
        5, 4, 3, 2,
        2, 2, 2, 2,
        3, 4, 5,
        5, 5, 5,
        4, 3,
        3, 3,
        4,
        4};
      int[] col = {0, 1, 2, 3, 4, 5, 6, 7,
        7, 7, 7, 7, 7, 7, 7,
        6, 5, 4, 3, 2, 1, 0,
        0, 0, 0, 0, 0, 0,
        1, 2, 3, 4, 5, 6,
        6, 6, 6, 6, 6,
        5, 4, 3, 2, 1,
        1, 1, 1, 1,
        2, 3, 4, 5,
        5, 5, 5,
        4, 3, 2,
        2, 2,
        3, 4,
        4,
        3};
      //lookUpSwap(row, col, blockSide);
      System.arraycopy(row, 0, refRow, 0, row.length);
      System.arraycopy(col, 0, refCol, 0, col.length);
    } else { //if (blockSide == 8) {
      int[] col = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
        14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
        13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
        13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
        12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
        12, 12, 12, 12, 12, 12, 12, 12, 12,
        11, 10, 9, 8, 7, 6, 5, 4, 3,
        3, 3, 3, 3, 3, 3, 3, 3,
        4, 5, 6, 7, 8, 9, 10, 11,
        11, 11, 11, 11, 11, 11, 11,
        10, 9, 8, 7, 6, 5, 4,
        4, 4, 4, 4, 4, 4,
        5, 6, 7, 8, 9, 10,
        10, 10, 10, 10, 10,
        9, 8, 7, 6, 5,
        5, 5, 5, 5,
        6, 7, 8, 9,
        9, 9, 9,
        8, 7, 6,
        6, 6,
        7, 8,
        8,
        7};

      int[] row = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
        14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
        13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
        13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
        12, 11, 10, 9, 8, 7, 6, 5, 4, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        4, 5, 6, 7, 8, 9, 10, 11, 12,
        12, 12, 12, 12, 12, 12, 12, 12, 12,
        11, 10, 9, 8, 7, 6, 5, 4,
        4, 4, 4, 4, 4, 4, 4, 4,
        5, 6, 7, 8, 9, 10, 11,
        11, 11, 11, 11, 11, 11, 11,
        10, 9, 8, 7, 6, 5,
        5, 5, 5, 5, 5, 5,
        6, 7, 8, 9, 10,
        10, 10, 10, 10, 10,
        9, 8, 7, 6,
        6, 6, 6, 6,
        7, 8, 9,
        9, 9, 9,
        8, 7,
        7, 7,
        8,
        8};
      //lookUpSwap(row, col, blockSide);
      System.arraycopy(row, 0, refRow, 0, row.length);
      System.arraycopy(col, 0, refCol, 0, col.length);
    }
  }

  float[] swapQuad(float[] array, int blockSide) {
    // Quadrants of the filter need to be swapped, as it's easier to program
    // filters to work on "edge" pixels representing high frequency data.
    int FFT_Side = 2 * blockSide;
    int tot = (int) Math.pow(FFT_Side, 2);
    int half_tot = tot / 2;
    float temp;

    String quad = "Top Left";
    int column = 1;
    // Swap quadrants
    for (int i = 0; i < half_tot; i++, column++) {

      // Operate on the quadrant variable to check where we are
      if (column > blockSide) {
        quad = "Top Right";
      }
      if (column > FFT_Side) {
        column = 1;
        quad = "Top Left";
      }

      // We should now have a variable that gives a top-right or -left
      //  answer
      int shift = 0;

      // Choose quadrant
      if (quad == "Top Left") {
        shift = i + half_tot + blockSide;
      } else if (quad == "Top Right") {
        shift = i + half_tot - blockSide;
      }

      // Assign new values
      temp = array[i];
      array[i] = array[shift];
      array[shift] = temp;
    }
    return array;

  }

  void customFilter(FHT fht, float[] FilterOriginal, int blockSide) {

    // Get pixels to operate with
    float[] fhtPixels = (float[]) fht.getPixels();
    // Create a new array for the filter to allow a quadrant swap without
    //  corrupting the original data
    float[] Filter = new float[FilterOriginal.length];
    System.arraycopy(FilterOriginal, 0, Filter, 0, Filter.length);
    // Swap the quandrants
    Filter = swapQuad(Filter, blockSide);

    for (int i = 0; i < fhtPixels.length; i++) {
      fhtPixels[i] = fhtPixels[i] * Filter[i];
    }
  }

  String progressReport() {
    // Create a tempory file to write data to in the temp dir
    String tempDir = "" + System.getProperty("java.io.tmpdir");
    tempDir += java.io.File.separator;
    tempDir += "KNM_prog.txt";

    // Initialise zero progress
    int progress = 0;

    canUpdate = true;
    updateProgress(tempDir, progress);

    IJ.doCommand("progress window");
    return tempDir;
  }

  void updateProgress(String fileName, int progress) {
    if (canUpdate) {
      FileOutputStream out; // declare a file output object
      PrintStream p; // declare a print stream object

      try {
        // Create a new file output stream
        out = new FileOutputStream(fileName);
        // Connect print stream to the output stream
        p = new PrintStream(out);
        p.println("FBNR progress...");
        p.println(Integer.toString(progress));
        p.close();
      } catch (Exception e) {
        IJ.showMessage("Error writing to file");
      }
      // Reset the ability to update
      canUpdate = false;
    }
  }
  Timer timer = new Timer(1000, new ActionListener() {

    public void actionPerformed(ActionEvent evt) {
      canUpdate = true;
    }
  });
}

class variance {
  // Object to hold variance data

  public double tot, tot0, res, signal, noise, max, scale;
  public boolean scaled;

  public variance() {
    // Initialise parameters to zero

    // Total variance of a block
    tot0 = 0;

    // Block variance post filtering
    tot = 0;

    // The difference between tot and tot2
    res = 0;

    // The noise variance, just equal to the mean pixel value
    noise = 0;

    // var.tot = var.signal - var.noise
    signal = 0;

    // The max pixel value of the unadulterated block
    max = 0;

    // Variable to hold whether the image has been scaled or not
    scaled = false;

    // How much the block has been scaled by
    scale = (double) 1;
  }
}
