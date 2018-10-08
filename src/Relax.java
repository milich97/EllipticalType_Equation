import java.util.ArrayList;

/**
 * Created by Миша on 01.10.2018.
 */
public class Relax extends Main {

    public Relax(double[][] Ustart) {
        int N = Ustart.length - 1;
        double xi = (Math.pow(Math.tan(Math.PI * h / 2), 2)) / 2;
        double ro = (1 - xi) / (1 + xi);
        int n = (int) (1 + Math.log(1 / eps) / Math.sqrt(xi));
        System.out.println();
        System.out.println("                                                      МЕТОД ВЕРХНЕЙ РЕЛАКСАЦИИ");
        System.out.println("Оценка числа итераций: " + n);
        System.out.println("Спектральный радиус матрицы перехода: " + ro);
        int k = 0;
        double omega = 2 / (1 + Math.sqrt(1 - Math.pow(ro, 2)));
        double[][] U1 = new double[Ustart.length][Ustart.length];
        for (int i = 0; i < Ustart.length; i++) {
            for (int j = 0; j < Ustart.length; j++) {
                U1[i][j] = Ustart[i][j];
            }
        }
        double[][] U0 = new double[Ustart.length][Ustart.length];
        for (int i = 0; i < Ustart.length; i++) {
            for (int j = 0; j < Ustart.length; j++) {
                U0[i][j] = Ustart[i][j];
            }
        }
        System.out.println("k      ||F-AU(k)||   rel.d.   ||U(k)-u*||   rel.error   ||U(k)-U(k-1)||   post.est.   app.spectral radius");

        ArrayList<Double> s = new ArrayList();
        double num, den;
        while ((norm2(U1, U_acc) / norm2(Ustart, U_acc)) > eps) {
            k++;


            for (int i = 1; i < U0.length - 1; i++) {
                for (int j = 1; j < U0.length - 1; j++) {
                    num = omega * (F[i][j] + 1 / h / h * (p.applyAsDouble(x[i] + h / 2, x[j]) * (U0[i + 1][j] - U0[i][j]) -
                            p.applyAsDouble(x[i] - h / 2, x[j]) * (U0[i][j] - U1[i - 1][j]) +
                            q.applyAsDouble(x[i], x[j] + h / 2) * (U0[i][j + 1] - U0[i][j]) -
                            q.applyAsDouble(x[i], x[j] - h / 2) * (U0[i][j] - U1[i][j - 1])));
                    den = 1 / h / h * (p.applyAsDouble(x[i] - h / 2, x[j]) +
                            p.applyAsDouble(x[i] + h / 2, x[j]) +
                            q.applyAsDouble(x[i], x[j] - h / 2) +
                            q.applyAsDouble(x[i], x[j] + h / 2));
                    U1[i][j] = U0[i][j] + num / den;
                }
            }

            s.add(norm2(U1, U0));
            if (k > 2) {
                System.out.println("        " + r(Math.sqrt(s.get(s.size() - 1) / s.get(s.size() - 3))));
            } else System.out.println();
            System.out.print(r1(k) + "        " + r(norm1(U1, F)) + "        " + r(norm1(U1, F) / norm1(Ustart, F)) +
                    "       " + r(norm2(U1, U_acc)) + "         " + r(norm2(U1, U_acc) / norm2(Ustart, U_acc)) +
                    "            " + r(norm2(U1, U0)) + "           " + r(ro * norm2(U1, U0) / (1 - ro)));

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
