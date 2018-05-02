package fastmath.java;

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
}
