
/**
 * Created by Миша on 08.10.2018.
 */
public class VariableTriangleMethod extends Main {
    public VariableTriangleMethod(double[][] Ustart) {
        int N = Ustart.length - 1;
        double xi = (Math.pow(Math.tan(Math.PI * h / 2), 2)) / 2;
        double ro = (1 - xi) / (1 + xi);
        int n = (int) (Math.log(1 / eps) / Math.log(1 / ro));
        double[][] U0 = new double[Ustart.length][Ustart.length];
        for (int i = 0; i < Ustart.length; i++) {
            for (int j = 0; j < Ustart.length; j++) {
                U0[i][j] = Ustart[i][j];
            }
        }
        double[][] U1 = new double[U0.length][U0.length];
        double[][] Phi = new double[N + 1][N + 1];
        double[][] W = new double[N + 1][N + 1];
        double[][] V = new double[N + 1][N + 1];


        System.out.println("                                                      ПОПЕРЕМЕННО-ТРЕУГОЛЬНЫЙ ИТЕРАЦИОННЫЙ МЕТОД");
        int k = 0;

        System.out.println("k       ||F-AU(k)||      rel.d.      ||U(k)-u*||     rel.error     ||U(k)-U(k-1)|| ");

        while ((norm2(U1, U_acc) / norm2(Ustart, U_acc)) > eps && k < n) {
            k++;
            double delta = 8 / h / h * Math.pow(Math.sin(Math.PI * h / 2), 2);
            double Delta = 16 / h / h;
            double eta = delta / Delta;
            double omega = 2 / Math.sqrt(delta * Delta);
            double gamma1 = delta / (2 + 2 * Math.sqrt(eta));
            double gamma2 = delta / (4 * Math.sqrt(eta));
            double tau = 2 / (gamma1 + gamma2);
            double kappa = omega / h / h;

            for (int i = 1; i < N; i++) {
                for (int j = 1; j < N; j++) {
                    Phi[i][j] = F[i][j] + 1 / h / h *
                            (p.applyAsDouble(x[i] + h / 2, x[j]) * (U0[i + 1][j] - U0[i][j]) -
                                    p.applyAsDouble(x[i] - h / 2, x[j]) * (U0[i][j] - U0[i - 1][j]) +
                                    q.applyAsDouble(x[i], x[j] + h / 2) * (U0[i][j + 1] - U0[i][j]) -
                                    q.applyAsDouble(x[i], x[j] - h / 2) * (U0[i][j] - U0[i][j - 1]));
                    W[i][j] = (kappa * (p.applyAsDouble(x[i] - h / 2, x[j]) * W[i - 1][j] +
                            q.applyAsDouble(x[i], x[j] - h / 2) * W[i][j - 1]) + Phi[i][j]) / (1 + 2 * kappa);
                }
            }
            for (int i = N - 1; i > 0; i--) {
                for (int j = N - 1; j > 0; j--) {
                    V[i][j] = (kappa * (p.applyAsDouble(x[i] + h / 2, x[j]) * V[i + 1][j] +
                            q.applyAsDouble(x[i], x[j] + h / 2) * V[i][j + 1]) + W[i][j]) / (1 + 2 * kappa);
                }
            }
            for (int i = 0; i <= N; i++) {
                for (int j = 0; j <= N; j++) {
                    U1[i][j] = U0[i][j] + tau * V[i][j];
                }
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


}
