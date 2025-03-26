package fastmath.java;

// https://github.com/lorentzenchr/scipy/blob/88abb5ed1cd85b23c4ae669d331a04599c28c145/scipy/optimize/_pava/pava_pybind.cpp

public final class Monotone {

    public static void pava(double[] x, double[] w, boolean asc) {
        boolean desc = !asc;
        int n = x.length;
        int[] r = new int[n+1];

        r[0]=0;
        r[1]=1;

        int b = 0;

        double xbp = x[0];
        double wbp = w[0];

        for(int i = 1; i < n; i++) {
            b++;

            double xb = x[i];
            double wb = w[i];
            double sb = 0;
            
            if ((asc && (xbp >= xb)) || (desc && (xbp <= xb))) {
                b--;

                sb = wbp * xbp + wb * xb;
                wb += wbp;
                xb = sb / wb;

                while (i<(n-1) && ((asc && (xb >= x[i+1])) || (desc && (xb <= x[i+1])))) {
                    i++;

                    sb += w[i] * x[i];
                    wb += w[i];
                    xb = sb / wb;
                }

                while (b>0 && ((asc && (x[b-1] >= xb)) || (desc && (x[b-1] <= xb)))) {
                    b--;

                    sb += w[b] * x[b];
                    wb += w[b];
                    xb = sb / wb;
                }
            }

            x[b] = xbp = xb;
            w[b] = wbp = wb;
            r[b+1] = i+1;
        }

        int f = n-1;

        for(int k = b; k >= 0; k--) {
            int t = r[k];
            double xk = x[k];

            for(int i = f; i >= t; i--) {
                x[i] = xk;
            }

            f = t - 1;
        }
        
    }
}
