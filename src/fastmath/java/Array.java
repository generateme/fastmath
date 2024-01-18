package fastmath.java;

import java.util.ArrayList;

public final class Array {

    public static int inc(int[] a, int idx) { return ++a[idx]; }
    public static long inc(long[] a, int idx) { return ++a[idx]; }
    public static double inc(double[] a, int idx) { return ++a[idx]; }

    public static int dec(int[] a, int idx) { return --a[idx]; }
    public static long dec(long[] a, int idx) { return --a[idx];}
    public static double dec(double[] a, int idx) { return --a[idx]; }

    public static int get(int[] a, int idx) { return a[idx]; }
    public static long get(long[] a, int idx) { return a[idx]; }
    public static double get(double[] a, int idx) { return a[idx]; }

    public static int set(int[] a, int idx, int v) { a[idx]=v; return v; }
    public static long set(long[] a, int idx, long v) { a[idx]=v; return v; }
    public static double set(double[] a, int idx, double v) { a[idx]=v; return v; }

    public static int add(int[] a, int idx, int v) { a[idx]+=v; return a[idx]; }
    public static long add(long[] a, int idx, long v) { a[idx]+=v; return a[idx]; }
    public static double add(double[] a, int idx, double v) { a[idx]+=v; return a[idx]; }

    public static int sub(int[] a, int idx, int v) { a[idx]-=v; return a[idx]; }
    public static long sub(long[] a, int idx, long v) { a[idx]-=v; return a[idx]; }
    public static double sub(double[] a, int idx, double v) { a[idx]-=v; return a[idx]; }

    public static int mult(int[] a, int idx, int v) { a[idx]*=v; return a[idx]; }
    public static long mult(long[] a, int idx, long v) { a[idx]*=v; return a[idx]; }
    public static double mult(double[] a, int idx, double v) { a[idx]*=v; return a[idx]; }

    public static int div(int[] a, int idx, int v) { a[idx]/=v; return a[idx]; }
    public static long div(long[] a, int idx, long v) { a[idx]/=v; return a[idx]; }
    public static double div(double[] a, int idx, double v) { a[idx]/=v; return a[idx]; }

    // 2d
    
    public static int inc2d(int[] a, int cols, int x, int y) { int idx=y*cols+x; return ++a[idx]; }
    public static long inc2d(long[] a, int cols, int x, int y) { int idx=y*cols+x; return ++a[idx]; }
    public static double inc2d(double[] a, int cols, int x, int y) { int idx=y*cols+x; return ++a[idx]; }

    public static int dec2d(int[] a, int cols, int x, int y) { int idx=y*cols+x; return --a[idx]; }
    public static long dec2d(long[] a, int cols, int x, int y) { int idx=y*cols+x; return --a[idx];}
    public static double dec2d(double[] a, int cols, int x, int y) { int idx=y*cols+x; return --a[idx]; }

    public static int get2d(int[] a, int cols, int x, int y) { int idx=y*cols+x; return a[idx]; }
    public static long get2d(long[] a, int cols, int x, int y) { int idx=y*cols+x; return a[idx]; }
    public static double get2d(double[] a, int cols, int x, int y) { int idx=y*cols+x; return a[idx]; }

    public static int set2d(int[] a, int cols, int x, int y, int v) { int idx=y*cols+x; a[idx]=v; return v; }
    public static long set2d(long[] a, int cols, int x, int y, long v) { int idx=y*cols+x; a[idx]=v; return v; }
    public static double set2d(double[] a, int cols, int x, int y, double v) { int idx=y*cols+x; a[idx]=v; return v; }

    public static int add2d(int[] a, int cols, int x, int y, int v) { int idx=y*cols+x; a[idx]+=v; return a[idx]; }
    public static long add2d(long[] a, int cols, int x, int y, long v) { int idx=y*cols+x; a[idx]+=v; return a[idx]; }
    public static double add2d(double[] a, int cols, int x, int y, double v) { int idx=y*cols+x; a[idx]+=v; return a[idx]; }

    public static int sub2d(int[] a, int cols, int x, int y, int v) { int idx=y*cols+x; a[idx]-=v; return a[idx]; }
    public static long sub2d(long[] a, int cols, int x, int y, long v) { int idx=y*cols+x; a[idx]-=v; return a[idx]; }
    public static double sub2d(double[] a, int cols, int x, int y, double v) { int idx=y*cols+x; a[idx]-=v; return a[idx]; }

    public static int mult2d(int[] a, int cols, int x, int y, int v) { int idx=y*cols+x; a[idx]*=v; return a[idx]; }
    public static long mult2d(long[] a, int cols, int x, int y, long v) { int idx=y*cols+x; a[idx]*=v; return a[idx]; }
    public static double mult2d(double[] a, int cols, int x, int y, double v) { int idx=y*cols+x; a[idx]*=v; return a[idx]; }

    public static int div2d(int[] a, int cols, int x, int y, int v) { int idx=y*cols+x; a[idx]/=v; return a[idx]; }
    public static long div2d(long[] a, int cols, int x, int y, long v) { int idx=y*cols+x; a[idx]/=v; return a[idx]; }
    public static double div2d(double[] a, int cols, int x, int y, double v) { int idx=y*cols+x; a[idx]/=v; return a[idx]; }

    public static double[][] mat2array2d(double a00, double a01, double a02, double a03,
                                         double a10, double a11, double a12, double a13,
                                         double a20, double a21, double a22, double a23,
                                         double a30, double a31, double a32, double a33) {
        return new double[][] {{a00,a01,a02,a03},
                               {a10,a11,a12,a13},
                               {a20,a21,a22,a23},
                               {a30,a31,a32,a33}};
    }

    public static double[][] mat2array2d(double a00, double a01, double a02,
                                         double a10, double a11, double a12,
                                         double a20, double a21, double a22) {
        return new double[][] {{a00,a01,a02},
                               {a10,a11,a12},
                               {a20,a21,a22}};
    }

    public static double[][] mat2array2d(double a00, double a01,
                                         double a10, double a11) {
        return new double[][] {{a00,a01},
                               {a10,a11}};
    }

    public static float[][] mat2array2d(float a00, float a01, float a02, float a03,
                                        float a10, float a11, float a12, float a13,
                                        float a20, float a21, float a22, float a23,
                                        float a30, float a31, float a32, float a33) {
        return new float[][] {{a00,a01,a02,a03},
                              {a10,a11,a12,a13},
                              {a20,a21,a22,a23},
                              {a30,a31,a32,a33}};
    }

    public static float[][] mat2array2d(float a00, float a01, float a02,
                                        float a10, float a11, float a12,
                                        float a20, float a21, float a22) {
        return new float[][] {{a00,a01,a02},
                              {a10,a11,a12},
                              {a20,a21,a22}};
    }

    public static float[][] mat2array2d(float a00, float a01,
                                        float a10, float a11) {
        return new float[][] {{a00,a01},
                              {a10,a11}};
    }

    public static double[] mat2array(double a00, double a01, double a02, double a03,
                                     double a10, double a11, double a12, double a13,
                                     double a20, double a21, double a22, double a23,
                                     double a30, double a31, double a32, double a33) {
        return new double[] {a00,a01,a02,a03,
                             a10,a11,a12,a13,
                             a20,a21,a22,a23,
                             a30,a31,a32,a33};
    }

    public static double[] mat2array(double a00, double a01, double a02,
                                     double a10, double a11, double a12,
                                     double a20, double a21, double a22) {
        return new double[] {a00,a01,a02,
                             a10,a11,a12,
                             a20,a21,a22};
    }

    public static double[] mat2array(double a00, double a01,
                                     double a10, double a11) {
        return new double[] {a00,a01,a10,a11};
    }

    public static float[] mat2array(float a00, float a01, float a02, float a03,
                                    float a10, float a11, float a12, float a13,
                                    float a20, float a21, float a22, float a23,
                                    float a30, float a31, float a32, float a33) {
        return new float[] {a00,a01,a02,a03,
                            a10,a11,a12,a13,
                            a20,a21,a22,a23,
                            a30,a31,a32,a33};
    }

    public static float[] mat2array(float a00, float a01, float a02,
                                    float a10, float a11, float a12,
                                    float a20, float a21, float a22) {
        return new float[] {a00,a01,a02,
                            a10,a11,a12,
                            a20,a21,a22};
    }

    public static float[] mat2array(float a00, float a01,
                                    float a10, float a11) {
        return new float[] {a00,a01,a10,a11};
    }

    public static ArrayList<double[]> mat2cols(double[][] m) {
        int rows = m.length;
        int cols = m[0].length;

        ArrayList<double[]> l = new ArrayList<double[]>();

        for(int c=0;c<cols;c++) {
            double[] t = new double[rows];
            for(int r=0;r<rows;r++) {
                t[r] = m[r][c];
            }
            l.add(t);
        }

        return l;
    }

    public static double[] mat2column(double[][] m, int col) {
        int rows = m.length;

        double[] t = new double[rows];
        for(int r=0;r<rows;r++) {
            t[r] = m[r][col];
        }

        return t;
    }

    public static double[] mat2diag(double[][] m) {
        int rows = m.length;

        double[] t = new double[rows];
        for(int r=0;r<rows;r++) {
            t[r] = m[r][r];
        }

        return t;
    }

    public static double[][] matadd(double[][] m1, double[][] m2) {
        int rows = m1.length;
        int cols = m1[0].length;

        double[][] t = new double[rows][cols];

        for(int c=0;c<cols;c++) {
            for(int r=0;r<rows;r++) {
                t[r][c] = m1[r][c] + m2[r][c];
            }
        }

        return t;
    }

    public static double[][] matadds(double[][] m, double s) {
        int rows = m.length;
        int cols = m[0].length;

        double[][] t = new double[rows][cols];

        for(int c=0;c<cols;c++) {
            for(int r=0;r<rows;r++) {
                t[r][c] = m[r][c] + s;
            }
        }

        return t;
    }
    
    public static double[][] matsub(double[][] m1, double[][] m2) {
        int rows = m1.length;
        int cols = m1[0].length;

        double[][] t = new double[rows][cols];

        for(int c=0;c<cols;c++) {
            for(int r=0;r<rows;r++) {
                t[r][c] = m1[r][c] - m2[r][c];
            }
        }

        return t;
    }

    public static double[][] matsub(double[][] m1) {
        int rows = m1.length;
        int cols = m1[0].length;

        double[][] t = new double[rows][cols];

        for(int c=0;c<cols;c++) {
            for(int r=0;r<rows;r++) {
                t[r][c] = -m1[r][c];
            }
        }

        return t;
    }

    public static double[][] matemulm(double[][] m1, double[][] m2) {
        int rows = m1.length;
        int cols = m1[0].length;

        double[][] t = new double[rows][cols];

        for(int c=0;c<cols;c++) {
            for(int r=0;r<rows;r++) {
                t[r][c] = m1[r][c] * m2[r][c];
            }
        }

        return t;
    }

    public static double[][] matmuls(double[][] m, double s) {
        int rows = m.length;
        int cols = m[0].length;

        double[][] t = new double[rows][cols];

        for(int c=0;c<cols;c++) {
            for(int r=0;r<rows;r++) {
                t[r][c] = m[r][c] * s;
            }
        }

        return t;
    }

}
