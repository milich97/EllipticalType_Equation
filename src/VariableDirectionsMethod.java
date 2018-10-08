
/**
 * Created by Миша on 08.10.2018.
 */
public class VariableDirectionsMethod extends Main {
    public VariableDirectionsMethod(double[][] Ustart) {
        int N = Ustart.length - 1;
        double[][] U0 = new double[Ustart.length][Ustart.length];
        for (int i = 0; i < Ustart.length; i++) {
            for (int j = 0; j < Ustart.length; j++) {
                U0[i][j] = Ustart[i][j];
            }
        }
        double[][] U1 = new double[U0.length][U0.length];
        double[][] Umid = new double[U0.length][U0.length];
        double[] A = new double[U0.length];
        double[] B = new double[U0.length];
        double[] C = new double[U0.length];
        double[] G = new double[U0.length];
        double[] help;


        System.out.println("                                                      МЕТОД ПЕРЕМЕННЫХ НАПРАВЛЕНИЙ");
        int k = 0;

        System.out.println("k       ||F-AU(k)||      rel.d.      ||U(k)-u*||     rel.error     ||U(k)-U(k-1)|| ");

        while ((norm2(U1, U_acc) / norm2(Ustart, U_acc)) > eps && k < 10) {
            k++;
            double delta = 4 / h / h * Math.pow(Math.sin(Math.PI * h / 2), 2);
            double Delta = 12 / h / h * Math.pow(Math.cos(Math.PI * h / 2), 2);
            double tau = 4 / Math.sqrt(delta * Delta);
            for (int i = 0; i < N + 1; i++) {
                Umid[i][0] = u_acc.applyAsDouble(x[i], 0);
            }
            for (int j = 1; j < N; j++) {
                for (int i = 1; i < N; i++) {
                    A[i] = tau / 2 / h / h * p.applyAsDouble(x[i] - h / 2, x[j]);
                    B[i] = tau / 2 / h / h * (p.applyAsDouble(x[i] - h / 2, x[j]) + p.applyAsDouble(x[i] + h / 2, x[j])) + 1;
                    C[i] = tau / 2 / h / h * p.applyAsDouble(x[i] + h / 2, x[j]);
                    G[i] = -U0[i][j] - tau / 2 * (f.applyAsDouble(x[i], x[j]) +
                            q.applyAsDouble(x[i], x[j] + h / 2) / h / h * (U0[i][j + 1] - U0[i][j]) -
                            q.applyAsDouble(x[i], x[j] - h / 2) / h / h * (U0[i][j] - U0[i][j - 1]));
                }
                A[0] = 0;
                A[N] = 0;
                B[0] = -1;
                B[N] = -1;
                C[0] = 0;
                C[N] = 0;
                G[0] = u_acc.applyAsDouble(0, x[j]);
                G[N] = u_acc.applyAsDouble(x[N], x[j]);

                help = strange_runningMethod1(A, B, C, G, j);
                for (int i = 0; i < N + 1; i++) {
                    Umid[i][j] = help[i];
                }
            }
            for (int i = 0; i < N + 1; i++) {
                Umid[i][N] = u_acc.applyAsDouble(x[i], x[N]);
            }

            for (int j = 0; j < N + 1; j++) {
                U1[0][j] = u_acc.applyAsDouble(0, x[j]);
            }
            for (int i = 1; i < N; i++) {
                for (int j = 1; j < N; j++) {
                    A[j] = tau / 2 / h / h * q.applyAsDouble(x[i], x[j] - h / 2);
                    B[j] = tau / 2 / h / h * (q.applyAsDouble(x[i], x[j] - h / 2) + q.applyAsDouble(x[i], x[j] + h / 2)) + 1;
                    C[j] = tau / 2 / h / h * q.applyAsDouble(x[i], x[j] + h / 2);
                    G[j] = -Umid[i][j] - tau / 2 * (f.applyAsDouble(x[i], x[j]) +
                            p.applyAsDouble(x[i] + h / 2, x[j]) / h / h * (Umid[i + 1][j] - Umid[i][j]) -
                            p.applyAsDouble(x[i] - h / 2, x[j]) / h / h * (Umid[i][j] - Umid[i - 1][j]));
                }
                A[0] = 0;
                A[N] = 0;
                B[0] = -1;
                B[N] = -1;
                C[0] = 0;
                C[N] = 0;
                G[0] = u_acc.applyAsDouble(x[i], 0);
                G[N] = u_acc.applyAsDouble(x[i], x[N]);

                help = strange_runningMethod2(A, B, C, G, i);
                for (int j = 0; j < N + 1; j++) {
                    U1[i][j] = help[j];
                }
            }
            for (int j = 0; j < N + 1; j++) {
                U1[N][j] = u_acc.applyAsDouble(x[N], x[j]);
            }


            System.out.println(r1(k) + "        " + r(norm1(U1, F)) + "        " + r(norm1(U1, F) / norm1(Ustart, F)) +
                    "       " + r(norm2(U1, U_acc)) + "         " + r(norm2(U1, U_acc) / norm2(Ustart, U_acc)) +
                    "            " + r(norm2(U1, U0)));

            for (int i = 0; i < U0.length; i++) {
                for (int j = 0; j < U0.length; j++) {
                    U0[i][j] = U1[i][j];
                }
            }
        }
        System.out.println();
        System.out.println();
        System.out.println("Значения U(k) на крупной сетке");
        for (int i = 0; i < U1.length; i += N / 5) {
            for (int j = 0; j < U1.length; j += N / 5) {
                System.out.print(r(U1[i][j]) + "    ");
            }
            System.out.println();
        }
        System.out.println();
        System.out.println("Значения u_acc на крупной сетке");
        for (int i = 0; i < U1.length; i += N / 5) {
            for (int j = 0; j < U1.length; j += N / 5) {
                System.out.print(r(u_acc.applyAsDouble(((double) i) / (U1.length - 1), ((double) j) / (U1.length - 1))) + "    ");

            }
            System.out.println();
        }


    }

    private double[] strange_runningMethod1(double[] A, double[] B, double[] C, double[] G, int j) {
        double[] y = new double[A.length];
        y[0] = u_acc.applyAsDouble(0, x[j]);
        y[A.length - 1] = u_acc.applyAsDouble(x[A.length - 1], x[j]);
        double[] help1;
        double[] a = new double[A.length - 2];
        double[] b = new double[A.length - 2];
        double[] c = new double[A.length - 2];
        double[] g = new double[A.length - 2];
        a[0] = 0;
        b[0] = B[1];
        c[0] = C[1];
        g[0] = G[1] - A[1] * u_acc.applyAsDouble(0, x[j]);
        for (int i = 0; i < A.length - 1 - 3; i++) {
            a[i + 1] = A[i + 2];
            b[i + 1] = B[i + 2];
            c[i + 1] = C[i + 2];
            g[i + 1] = G[i + 2];
        }
        a[A.length - 3] = A[A.length - 1 - 1];
        b[A.length - 3] = B[A.length - 2];
        c[A.length - 3] = 0;
        g[A.length - 3] = G[A.length - 2] - C[A.length - 2] * u_acc.applyAsDouble(x[A.length - 1], x[j]);

        help1 = runningMethod(a, b, c, g);
        for (int i = 0; i < A.length - 1 - 1; i++) {
            y[i + 1] = help1[i];
        }
        return y;
    }

    private double[] strange_runningMethod2(double[] A, double[] B, double[] C, double[] G, int i) {
        double[] y = new double[A.length];
        y[0] = u_acc.applyAsDouble(x[i], 0);
        y[A.length - 1] = u_acc.applyAsDouble(x[i], x[A.length - 1]);
        double[] help1;
        double[] a = new double[A.length - 2];
        double[] b = new double[A.length - 2];
        double[] c = new double[A.length - 2];
        double[] g = new double[A.length - 2];
        a[0] = 0;
        b[0] = B[1];
        c[0] = C[1];
        g[0] = G[1] - A[1] * u_acc.applyAsDouble(x[i], 0);
        for (int j = 0; j < A.length - 1 - 3; j++) {
            a[j + 1] = A[j + 2];
            b[j + 1] = B[j + 2];
            c[j + 1] = C[j + 2];
            g[j + 1] = G[j + 2];
        }
        a[A.length - 3] = A[A.length - 1 - 1];
        b[A.length - 3] = B[A.length - 2];
        c[A.length - 3] = 0;
        g[A.length - 3] = G[A.length - 2] - C[A.length - 2] * u_acc.applyAsDouble(x[i], x[A.length - 1]);

        help1 = runningMethod(a, b, c, g);
        for (int j = 0; j < A.length - 1 - 1; j++) {
            y[j + 1] = help1[j];
        }
        return y;
    }

    private static double[] runningMethod(double[] A, double[] B, double[] C, double[] G) {

/*
Имеем линейную замкнутую систему (n + 1)-го порядка:
- B_{0}*y_{0} + C_{0}*y_{1} = G_{0}
A_{i}*y_{i-1} - B_{i}*y_{i} + C_{i}*y_{i+1} = G_{i}
. . . . . . . . . . . . . . . . . . . . . . .
A_{n}*y_{n-1} - B_{n}*y_{n} = G_{n}
*/
        // A - массив коэфиниентов перед y_{i-1}
        // B - массив коэфиниентов перед (-y_{i})
        // C - массив коэфиниентов перед y_{i+1}
        // G - массив свободных членов

/*
Ищем решение в виде:
y_{i} = s_{i} * y_{i+1} + t_{i}, i = 0,...,n-1
*/
        double[] s = new double[A.length]; // массив значений s_{i}
        double[] t = new double[G.length]; // массив значений t_{i}
        double[] y = new double[A.length]; // массив значений в узлах сетки


        double s_0 = C[0] / B[0];
        double t_0 = -G[0] / B[0];

        s[0] = s_0;
        for (int i = 1; i < A.length - 1; i++) {
            s[i] = C[i] / (B[i] - A[i] * s[i - 1]);
        }
        s[A.length - 1] = 0;

        t[0] = t_0;
        for (int i = 1; i < A.length; i++) {
            t[i] = (A[i] * t[i - 1] - G[i]) / (B[i] - A[i] * s[i - 1]);
        }

        y[A.length - 1] = t[A.length - 1];
        for (int i = A.length - 1 - 1; i > -1; i--) {
            y[i] = s[i] * y[i + 1] + t[i];
        }

        return y;
    }
}
