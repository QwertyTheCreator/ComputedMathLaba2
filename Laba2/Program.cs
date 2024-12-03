using System;

public class Program
{
    private const double EPSILON = 1e-5;
    private const int N = 15;
    private static readonly double[] X = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };

    // Метод для вычисления спектральной нормы матрицы
    public static double CalculateSpectralNorm(double[,] matrix)
    {
        double[] x = InitializeX(1);
        double[] y = new double[N];
        double lambdaPrev = 0;
        double lambda;
        double norm;

        do
        {
            for (int i = 0; i < N; i++)
            {
                y[i] = 0;
                for (int j = 0; j < N; j++)
                {
                    y[i] += matrix[i, j] * x[j];
                }
            }

            lambda = 0;
            for (int i = 0; i < N; i++)
            {
                lambda = Math.Max(lambda, Math.Abs(y[i]));
            }

            for (int i = 0; i < N; i++)
            {
                x[i] = y[i] / lambda;
            }

            norm = Math.Abs(lambda - lambdaPrev);
            lambdaPrev = lambda;
        } while (norm > EPSILON);

        return lambda;
    }

    // Вычисление нормы разности между точным результатом и найденным
    private static double FindNormResult(double[] x)
    {
        double normMax = 0.0;
        for (int i = 0; i < N; i++)
        {
            if (normMax < Math.Abs(X[i] - x[i]))
            {
                normMax = Math.Abs(X[i] - x[i]);
            }
        }
        return normMax;
    }

    // Инициализация x
    private static double[] InitializeX(double value)
    {
        double[] x = new double[N];
        for (int i = 0; i < N; i++)
        {
            x[i] = value;
        }
        return x;
    }

    // Решение системы методом простой итерации
    public static int SolveBySimpleIteration(double[,] B, double[] b, double[] x, double normB)
    {
        double[] xNew = new double[N];
        int numberIterations = 0;
        double normDiff;

        do
        {
            normDiff = 0.0;
            for (int i = 0; i < N; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < N; j++)
                {
                    sum += B[i, j] * x[j];
                }
                xNew[i] = b[i] + sum;
                normDiff += Math.Abs(xNew[i] - x[i]);
            }

            Array.Copy(xNew, x, N);
            numberIterations++;
        } while (normDiff * normB / (1 - normB) > EPSILON);

        Console.WriteLine("Метод простой итерации");
        Console.WriteLine("Решение найдено за " + numberIterations + " итераций.");
        Console.WriteLine("Решение:");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine($"x[{i}] = {x[i]}");
        }

        double normRes = FindNormResult(x);
        Console.WriteLine(normRes);

        return numberIterations;
    }

    // Решение системы методом Зейделя
    private static void SolveBySeidel(double[,] B, double[] b, double[] x, int numberIterations)
    {
        double[] xNew = new double[N];
        for (int k = 0; k < numberIterations; k++)
        {
            for (int i = 0; i < N; i++)
            {
                double sum = 0.0;
                for (int j = i; j < N; j++)
                {
                    sum += B[i, j] * x[j];
                }

                for (int j = 0; j < i; j++)
                {
                    sum += B[i, j] * xNew[j];
                }
                xNew[i] = b[i] + sum;
            }

            Array.Copy(xNew, x, N);
        }

        Console.WriteLine("Метод Зейделя");
        Console.WriteLine("Решение найдено за " + numberIterations + " итераций.");
        Console.WriteLine("Решение:");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine($"x[{i}] = {x[i]}");
        }

        double normRes = FindNormResult(x);
        Console.WriteLine(normRes);
    }

    // Решение системы методом минимальных невязок
    static (double[], int) MinResidualMethod(double[] x0, double[,] A, double[] b, int maxIterations = 1000)
    {
        double[] x = (double[])x0.Clone();
        double[] r = Multiply(A, x).Zip(b, (a, bVal) => a - bVal).ToArray(); //Вектор невязок

        for (int iterations = 0; iterations < maxIterations; iterations++)
        {
            double tao = Multiply(A, r).Zip(r, (a, bVal) => a * bVal).Sum() / Math.Pow(Norm(Multiply(A, r)), 2);
            double[] xNew = x.Zip(r, (xVal, rVal) => xVal - tao * rVal).ToArray();

            r = Multiply(A, xNew).Zip(b, (a, bVal) => a - bVal).ToArray();
            if (Norm(r) < EPSILON)
            {
                Console.WriteLine("Метод минимальных невязок");
                Console.WriteLine("Решение найдено за " + iterations + " итераций.");
                Console.WriteLine("Решение:");
                for (int i = 0; i < N; i++)
                {
                    Console.WriteLine($"x[{i}] = {x[i]}");
                }

                double normRes = FindNormResult(x);
                Console.WriteLine(normRes);

                return (xNew, iterations + 1);
            }
            x = xNew;
        }
        throw new InvalidOperationException("Метод не сошелся за заданное количество итераций");
    }

    static double Norm(double[] vector)
    {
        return Math.Sqrt(vector.Sum(v => v * v));
    }

    static double[] Multiply(double[,] A, double[] x)
    {
        double[] result = new double[N];
        for (int i = 0; i < N; i++)
        {
            result[i] = Enumerable.Range(0, N).Sum(j => A[i, j] * x[j]);
        }
        return result;
    }

    // Решение системы
    public static void Solve(double[,] A, double[] b, double[] x)
    {
        double[] bNew = new double[N];

        // Вычисление нормы матрицы A
        double normA = CalculateSpectralNorm(A);

        double e = 1;
        double normB = 0;
        double[,] B = new double[N, N];

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                B[i, j] = (i != j) ? -A[i, j] / normA
                                   : e - A[i, j] / normA;
            }

            bNew[i] = b[i] / normA;

            if (normB < Math.Abs(bNew[i]))
            {
                normB = Math.Abs(bNew[i]);
            }
        }

        // Вычисление нормы матрицы B
        double normBValue = CalculateSpectralNorm(B);

        double aprioriEstimate = Math.Log(((1 - normBValue) * EPSILON) / normB) / Math.Log(normBValue);
        int iterationsNeeded = (int)Math.Ceiling(aprioriEstimate);
        Console.WriteLine("Рассчитанное число итераций: " + iterationsNeeded);

        int numberIteration = SolveBySimpleIteration(B, bNew, x, normBValue);

        x = InitializeX(0);
        SolveBySeidel(B, bNew, x, numberIteration);

        x = InitializeX(0);
        MinResidualMethod(x, A, b); // Здесь заменен метод градиентного спуска
    }

    public static void Main(string[] args)
    {
        double[,] A = new double[N, N];
        double[] b = new double[N];
        double[] x = new double[N];

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (i == j)
                {
                    A[i, j] = 5 * (i + 1);
                    continue;
                }

                A[i, j] = 0.1 * (i + 1) * (j + 1);
            }
        }

        for (int i = 0; i < N; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < N; j++)
            {
                sum += A[i, j] * X[j];
            }
            b[i] = sum;
            x[i] = 0;
        }

        Solve(A, b, x);
    }
}