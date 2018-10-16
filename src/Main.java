import java.util.function.DoubleBinaryOperator;

/**
 * Created by Миша on 08.09.2018.
 */
public class Main {
    protected static DoubleBinaryOperator p = (x, y) -> 1 + 2 * x;
    protected static DoubleBinaryOperator q = (x, y) -> 1;
    protected static DoubleBinaryOperator u_acc = (x, y) -> x * x * y * y * (1 + y);
    protected static DoubleBinaryOperator f = (x, y) -> -(2 * (y * y + y * y * y) * (1 + 4 * x) + 2 * x * x * (3 * y + 1));
    protected static double[] x;
    protected static double h;
    protected static double eps = 0.005;
    protected static double[][] U_acc;
    protected static double[][] F;

    public static void main(String[] args) {
        int N = 10;

        h = 1 / ((double) N);

        x = new double[N + 1];
        for (int i = 0; i < x.length; i++) {
            x[i] = i * h;
        }

        double[][] U0 = new double[N + 1][N + 1];
        for (int i = 0; i < U0.length; i++) {
            U0[0][i] = u_acc.applyAsDouble(0, x[i]);
            U0[N][i] = u_acc.applyAsDouble(1, x[i]);
        }
        for (int i = 0; i < U0.length - 2; i++) {
            U0[i + 1][0] = u_acc.applyAsDouble(x[i + 1], 0);
            U0[i + 1][N] = u_acc.applyAsDouble(x[i + 1], 1);
        }

        F = new double[N + 1][N + 1];
        for (int i = 0; i < N + 1; i++) {
            for (int j = 0; j < N + 1; j++) {
                F[i][j] = f.applyAsDouble(x[i], x[j]);
            }
        }

        U_acc = new double[N + 1][N + 1];
        for (int i = 0; i < U_acc.length; i++) {
            for (int j = 0; j < U_acc.length; j++) {
                U_acc[i][j] = u_acc.applyAsDouble(x[i], x[j]);
            }
        }


        System.out.println("Мера аппроксимации ДУ на точном решении: " + r(norm1(U_acc, F)));
        System.out.println("Норма невязки нулевого приближения: " + r(norm1(U0, F)));
//        new IterationMethodModified(U0);
//        new Seidel(U0);
//        new Relax(U0);
        new VariableTriangleMethod(U0);
//        new VariableDirectionsMethod(U0);

    }

    protected static double norm1(double[][] u, double[][] F) {
        double max = 0;
        double value;
        for (int i = 1; i < F.length - 1; i++) {
            for (int j = 1; j < F.length - 1; j++) {
                value = Math.abs(F[i][j] +
                        p.applyAsDouble(x[i] + h / 2, x[j]) * (u[i + 1][j] - u[i][j]) / h / h -
                        p.applyAsDouble(x[i] - h / 2, x[j]) * (u[i][j] - u[i - 1][j]) / h / h +
                        q.applyAsDouble(x[i], x[j] + h / 2) * (u[i][j + 1] - u[i][j]) / h / h -
                        q.applyAsDouble(x[i], x[j] - h / 2) * (u[i][j] - u[i][j - 1]) / h / h);
                if (value > max) max = value;
            }
        }
        return max;
    }

    protected static double norm2(double[][] u, double[][] U) {
        double max = 0;
        double value;
        for (int i = 1; i < U.length - 1; i++) {
            for (int j = 1; j < U.length - 1; j++) {
                value = Math.abs(u[i][j] - U[i][j]);
                if (value > max) max = value;
            }
        }
        return max;
    }

    protected static String r(double value) {
        return String.format("%5.4f", value);
    }

    protected static String r1(int value) {
        if (value < 10) return value + " ";
        else return value + "";
    }
}
