import java.util.ArrayList;

/**
 * Created by Миша on 08.10.2018.
 */
public class VariableTriangleMethod extends Main {
    public VariableTriangleMethod(double[][] Ustart) {
        int N = Ustart.length - 1;
        double xi = (Math.pow(Math.tan(Math.PI * h / 2), 2)) / 2;
        double ro = (1 - xi) / (1 + xi);
        int n = (int) (Math.log(1 / eps) / Math.log(1 / ro));
        double[][] u0 = new double[Ustart.length][Ustart.length];
        for (int i = 0; i < Ustart.length; i++) {
            for (int j = 0; j < Ustart.length; j++) {
                u0[i][j] = Ustart[i][j];
            }
        }
        double[][] u1 = new double[u0.length][u0.length];
        double[][] Phi = new double[N + 1][N + 1];
        double[][] U_ = new double[N + 1][N + 1];
        double[][] U = new double[N + 1][N + 1];
        double delta = 8 / h / h * Math.pow(Math.sin(Math.PI * h / 2), 2);
        double Delta = 16 / h / h;
        double eta = delta / Delta;
        double omega = 2 / Math.sqrt(delta * Delta);
        double gamma1 = delta / (2 + 2 * Math.sqrt(eta));
        double gamma2 = delta / (4 * Math.sqrt(eta));
        double tau = 2 / (gamma1 + gamma2);
        double kappa = omega / h / h;

        System.out.println("                                                      ПОПЕРЕМЕННО-ТРЕУГОЛЬНЫЙ ИТЕРАЦИОННЫЙ МЕТОД");
        int k = 0;

        System.out.println("k       ||F-AU(k)||     rel.d.      ||U(k)-u*||    rel.error      ||U(k)-U(k-1)||    ");


        while ((norm2(u1, U_acc) / norm2(Ustart, U_acc)) > eps && k < 10) {
            k++;

            for (int i = 1; i < N; i++) {
                for (int j = 1; j < N; j++) {
                    Phi[i][j] = F[i][j] + 1 / h / h *
                            (p.applyAsDouble(x[i] + h / 2, x[j]) * (u0[i + 1][j] - u0[i][j]) -
                                    p.applyAsDouble(x[i] - h / 2, x[j]) * (u0[i][j] - u0[i - 1][j]) +
                                    q.applyAsDouble(x[i], x[j] + h / 2) * (u0[i][j + 1] - u0[i][j]) -
                                    q.applyAsDouble(x[i], x[j] - h / 2) * (u0[i][j] - u0[i][j - 1]));
                    U_[i][j] = (kappa * (p.applyAsDouble(x[i] - h / 2, x[j]) * U_[i - 1][j] +
                            q.applyAsDouble(x[i], x[j] - h / 2) * U_[i][j - 1]) + Phi[i][j]) / (1 +
                            p.applyAsDouble(x[i] - h / 2, x[j]) * kappa +
                            q.applyAsDouble(x[i], x[j] - h / 2) * kappa);

                }
            }
            for (int i = N - 1; i > 0; i--) {
                for (int j = N - 1; j > 0; j--) {
                    U[i][j] = (kappa * (p.applyAsDouble(x[i] + h / 2, x[j]) * U[i + 1][j] +
                            q.applyAsDouble(x[i], x[j] + h / 2) * U[i][j + 1]) + U_[i][j]) / (1 +
                            p.applyAsDouble(x[i] + h / 2, x[j]) * kappa +
                            q.applyAsDouble(x[i], x[j] + h / 2) * kappa);
                }
            }
            for (int i = 0; i <= N; i++) {
                for (int j = 0; j <= N; j++) {
                    u1[i][j] = u0[i][j] + tau * U[i][j];
                }
            }


            System.out.println(r1(k) + "        " + r(norm1(u1, F)) + "        " + r(norm1(u1, F) / norm1(Ustart, F)) +
                    "       " + r(norm2(u1, U_acc)) + "         " + r(norm2(u1, U_acc) / norm2(Ustart, U_acc)) +
                    "            " + r(norm2(u1, u0)));

            for (int i = 0; i < u0.length; i++) {
                for (int j = 0; j < u0.length; j++) {
                    u0[i][j] = u1[i][j];
                }
            }
        }
        System.out.println();
        System.out.println();
        System.out.println("Значения U(k) на крупной сетке");
        for (int i = 0; i < u1.length; i += N / 5) {
            for (int j = 0; j < u1.length; j += N / 5) {
                System.out.print(r(u1[i][j]) + "    ");
            }
            System.out.println();
        }
        System.out.println();
        System.out.println("Значения u_acc на крупной сетке");
        for (int i = 0; i < u1.length; i += N / 5) {
            for (int j = 0; j < u1.length; j += N / 5) {
                System.out.print(r(u_acc.applyAsDouble(((double) i) / (u1.length - 1), ((double) j) / (u1.length - 1))) + "    ");

            }
            System.out.println();
        }


    }


}
//                    W[i][j]=omega/h/h*(p.applyAsDouble(x[i]-h/2,x[j])*W[i-1][j]+q.applyAsDouble(x[i],x[j]-h/2)*W[i][j-1])+
//                            Phi[i][j]/(1+omega/h/h*(p.applyAsDouble(x[i]-h/2,x[j])+q.applyAsDouble(x[i],x[j]-h/2)));

//                V[i][j]=omega/h/h*(p.applyAsDouble(x[i]+h/2,x[j])*V[i+1][j]+q.applyAsDouble(x[i],x[j]+h/2)*V[i][j+1])+
//                        W[i][j]/(1+omega/h/h*(p.applyAsDouble(x[i]+h/2,x[j])+q.applyAsDouble(x[i],x[j]+h/2)));