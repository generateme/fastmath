package fastmath.java;

// https://github.com/lorentzenchr/scipy/blob/88abb5ed1cd85b23c4ae669d331a04599c28c145/scipy/optimize/_pava/pava_pybind.cpp

public final class Monotone {

    public static int pava_step1(double[] y, double[] w, int[] r, int order) {
        return pava_step1(null, y, w, r, order);
    }

    // order = 0: weak asc
    // order = 1: weak desc
    // order = 2: strict asc
    // order = 3: string desc
    
    public static int pava_step1(double[] x, double[] y, double[] w, int[] r, int order) {
        boolean desc = (order == 1 || order == 3);
        boolean asc = !desc;
        boolean strict = (order == 2 || order == 3);
        boolean weak = !strict;
        boolean do_x = (x != null);

        int n = y.length;

        r[0] = 0;
        r[1] = 1;

        int b = 0;

        double xbp = do_x ? x[0] : 0.0;
        double ybp = y[0];
        double wbp = w[0];

        for(int i = 1; i < n; i++) {
            b++;

            double xb = do_x ? x[i] : 0.0;
            double yb = y[i];
            double wb = w[i];
            double sby = 0.0;
            double sbx = 0.0;
            
            if ((weak && ((asc && (ybp >= yb)) || (desc && (ybp <= yb)))) ||
                (strict && ((asc && (ybp > yb)) || (desc && (ybp < yb))))) {
                
                b--;

                if (do_x) sbx = wbp * xbp + wb * xb;
                sby = wbp * ybp + wb * yb;
                 
                wb += wbp;

                if (do_x) xb = sbx / wb;
                yb = sby / wb;
                 
                while (i<(n-1) && ((weak && ((asc && (yb >= y[i+1])) || (desc && (yb <= y[i+1])))) ||
                                   (strict && ((asc && (yb > y[i+1])) || (desc && (yb < y[i+1])))))) {
                    
                    i++;

                    if (do_x) sbx += w[i] * x[i];
                    sby += w[i] * y[i];

                    wb += w[i];

                    if (do_x) xb = sbx / wb;
                    yb = sby / wb;
                }

                while (b>0 && ((weak && ((asc && (y[b-1] >= yb)) || (desc && (y[b-1] <= yb)))) ||
                               (strict && ((asc && (y[b-1] > yb)) || (desc && (y[b-1] < yb)))))) {
                    b--;

                    if (do_x) sbx += w[b] * x[b];
                    sby += w[b] * y[b];
                     
                    wb += w[b];

                    if (do_x) xb = sbx / wb;
                    yb = sby / wb;                 
                }
            }

            if (do_x) x[b] = xbp = xb;             
            y[b] = ybp = yb;
            w[b] = wbp = wb;
            r[b+1] = i+1;
        }

        return b;
    }

    public static void pava_step2(double[] y, int[] r, int b) {
        int f = y.length - 1;

        for(int k = b; k >= 0; k--) {
            int t = r[k];
            double yk = y[k];

            for(int i = f; i >= t; i--) {
                y[i] = yk;
            }

            f = t - 1;
        }   
    }
}
