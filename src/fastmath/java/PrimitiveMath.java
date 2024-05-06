// Primitive Math operations with all variants
package fastmath.java;

public final class PrimitiveMath {
    public static long add(long a, long b) { return a+b; }
    public static double add(double a, long b) { return a+b; }
    public static double add(long a, double b) { return a+b; }
    public static double add(double a, double b) { return a+b; }

    public static long subtract(long a, long b) { return a-b; }
    public static double subtract(double a, long b) { return a-b; }
    public static double subtract(long a, double b) { return a-b; }
    public static double subtract(double a, double b) { return a-b; }
    
    public static long negate(long a) { return -a; }
    public static double negate(double a) { return -a; }
    
    public static long multiply(long a, long b) { return a*b; }
    public static double multiply(double a, long b) { return a*b; }
    public static double multiply(long a, double b) { return a*b; }
    public static double multiply(double a, double b) { return a*b; }

    public static long divide(long a, long b) { return a/b; }
    public static double divide(double a, long b) { return a/b; }
    public static double divide(long a, double b) { return a/b; }
    public static double divide(double a, double b) { return a/b; }
    
    public static long shiftLeft(long a, long n) { return a << n; }
    public static long shiftRight(long a, long n) { return a >> n; }
    public static long unsignedShiftRight(long a, long n) { return a >>> n; }
    public static int unsignedShiftRight(int a, long n) { return a >>> n; }

    public static long bitAnd(long a, long b) { return a & b; }
    public static long bitNand(long a, long b) { return ~(a & b); }
    public static long bitAndNot(long a, long b) { return a & ~b; }
    public static long bitOr(long a, long b) { return a | b; }
    public static long bitNor(long a, long b) { return ~(a | b); }
    public static long bitXor(long a, long b) { return a ^ b; }
    public static long bitXNor(long a, long b) { return ~(a ^ b); }
    public static long bitNot(long a) { return ~a; }

    public static long bitSet (long a, long n) {return a | (1L << n);}
    public static long bitClear (long a, long n) {return a & ~(1L << n);}
    public static long bitFlip (long a, long n) {return a ^ (1L << n);}
    public static boolean bitTest (long a, long n) {return (a & (1L << n)) != 0L;}
    
    public static long inc(long a) { return a+1L; }
    public static double inc(double a) { return a+1.0; }
    public static long dec(long a) { return a-1L; }
    public static double dec(double a) { return a-1.0; }

    public static double reciprocal(long a) { return 1.0/a; }
    public static double reciprocal(double a) { return 1.0/a; }
    
    public static long remainder(long a, long b) { return a%b; }
    public static double remainder(double a, long b) { return a%b; }
    public static double remainder(long a, double b) { return a%b; }
    public static double remainder(double a, double b) { return a%b; }

    public static long modulus(long a, long b) { long t=a%b; return ((t==0L) || (a>0L)==(b>0L))?t:t+b; }
    public static double modulus(double a, long b) { double t=a%b; return ((t==0.0) || (a>0.0)==(b>0L))?t:t+b; }
    public static double modulus(long a, double b) { double t=a%b; return ((t==0.0) || (a>0L)==(b>0.0))?t:t+b; }
    public static double modulus(double a, double b) { double t=a%b; return ((t==0.0) || (a>0.0)==(b>0.0))?t:t+b; }

    public static long quotient(long a, long b) { return a/b; }
    public static double quotient(double a, long b) { return (double)((long)(a/b)); }
    public static double quotient(long a, double b) { return (double)((long)(a/b)); }
    public static double quotient(double a, double b) { return (double)((long)(a/b)); }

    public static boolean lt(long a, long b) { return a<b; }
    public static boolean lt(double a, long b) { return a<b; }
    public static boolean lt(long a, double b) { return a<b; }
    public static boolean lt(double a, double b) { return a<b; }

    public static boolean gt(long a, long b) { return a>b; }
    public static boolean gt(double a, long b) { return a>b; }
    public static boolean gt(long a, double b) { return a>b; }
    public static boolean gt(double a, double b) { return a>b; }

    public static boolean lte(long a, long b) { return a<=b; }
    public static boolean lte(double a, long b) { return a<=b; }
    public static boolean lte(long a, double b) { return a<=b; }
    public static boolean lte(double a, double b) { return a<=b; }

    public static boolean gte(long a, long b) { return a>=b; }
    public static boolean gte(double a, long b) { return a>=b; }
    public static boolean gte(long a, double b) { return a>=b; }
    public static boolean gte(double a, double b) { return a>=b; }

    public static boolean eq(long a) { return true; }
    public static boolean eq(double a) { return true; }

    public static boolean eq(long a, long b) { return a==b; }
    public static boolean eq(double a, long b) { return a==b; }
    public static boolean eq(long a, double b) { return a==b; }
    public static boolean eq(double a, double b) { return a==b; }

    public static boolean neq(long a, long b) { return a!=b; }
    public static boolean neq(double a, long b) { return a!=b; }
    public static boolean neq(long a, double b) { return a!=b; }
    public static boolean neq(double a, double b) { return a!=b; }

    public static boolean isZero(long a) { return a==0L; }
    public static boolean isZero(double a) { return a==0.0; }
    public static boolean isOne(long a) { return a==1L; }
    public static boolean isOne(double a) { return a==1.0; }

    public static boolean isNeg(long a) { return a<0L; }
    public static boolean isNeg(double a) { return a<0.0; }
    public static boolean isPos(long a) { return a>0L; }
    public static boolean isPos(double a) { return a>0.0; }

    public static boolean isNNeg(long a) { return a>=0L; }
    public static boolean isNNeg(double a) { return a>=0.0; }
    public static boolean isNPos(long a) { return a<=0L; }
    public static boolean isNPos(double a) { return a<=0.0; }

    public static boolean and(boolean a, boolean b) { return a && b; }
    public static boolean or(boolean a, boolean b) { return a || b; }
    public static boolean not(boolean a) { return !a; }
    public static boolean xor(boolean a, boolean b) { return (a || b) && !(a && b); }

    public static long min(long a, long b) { return a<b?a:b; }
    public static double min(double a, long b) { return a<b?a:b; }
    public static double min(long a, double b) { return a<b?a:b; }
    public static double min(double a, double b) { return a<b?a:b; }
    
    public static long max(long a, long b) { return a>b?a:b; }
    public static double max(double a, long b) { return a>b?a:b; }
    public static double max(long a, double b) { return a>b?a:b; }
    public static double max(double a, double b) { return a>b?a:b; }

    public static boolean isEven(long a) { return (a&1)==0; }
    public static boolean isOdd(long a) { return (a&1)==1; }

    public static double norm(double v, double start, double stop) {
        return start == stop ? (v <= start ? 0.0 : 1.0) : (v - start) / (stop - start);
    }
    public static double norm(double v, double start1, double stop1, double start2, double stop2) {
        return start2 + (stop2 - start2) * norm(v, start1, stop1);
    }

    public static double lerp(double a, double b, double t) { return a + t * (b - a); }
    public static double hermite(double v) { return v * v * (3.0 - 2.0 * v); }
    public static double quintic(double v) { return v * v * v * (v * (v * 6.0 - 15.0) + 10.0); }

    public static long fastFloor(double v) { return v > 0.0 ? (long)v : (long)v - 1L; }
    public static long fastCeil(double v) { return - fastFloor(-v); }
    public static long fastRound(double v) { return v > 0.0 ? (long)(v + 0.5) : (long)(v - 0.5); }
}
