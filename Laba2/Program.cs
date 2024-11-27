const int SIZE = 15;
const int GROUP_NUMBER = 2;
const double EPSILON = 1e-5;

var matrix = BuildMatrix(SIZE);
var solution = BuildStarterSolution(GROUP_NUMBER, SIZE);
var rightPart = RandomRightPart(SIZE);

var (B, d) = ConvertToIterativeForm(matrix, rightPart);

var 





static double[] MultiplyMatrixByVector(double[][] matrix, double[] vector)
{
    int rows = SIZE;
    int cols = SIZE;

    double[] result = new double[rows];

    for (int i = 0; i < rows; i++)
    {
        result[i] = 0;
        for (int j = 0; j < cols; j++)
        {
            result[i] += matrix[i][j] * vector[j];
        }
    }

    return result;
}

static (double[][], double[]) ConvertToIterativeForm(double[][] A, double[] b)
{
    int n = SIZE;
    double[][] B = new double[n][];
    double[] d = new double[n];

    for (int i = 0; i < n; i++)
    {
        B[i] = new double[n];
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                B[i][j] = -A[i][j] / A[i][i];
            }
        }

        d[i] = b[i] / A[i][i];
    }

    return (B, d);
}

static double[] BuildStarterSolution(int groupNumber = 2, int size = 15)
{
    var x = new double[SIZE];
    for(int i = 0; i < SIZE; i++)
    {
        x[i] = groupNumber + i;
    }

    return x;
}

static double[] RandomRightPart(int size = 15)
{
    var b = new double[SIZE];
    for (int i = 0; i < SIZE; i++)
    {
        b[i] = Random.Shared.NextDouble()*20;
    }

    return b;
}

static double[][] BuildMatrix(int size = 15)
{
    double[][] matrix = new double[size][];

    for (int i = 0; i < size; i++)
    {
        matrix[i] = new double[size];
        for (int j = 0; j < size; j++)
        {
            if (i == j)
            {
                matrix[i][j] = 5 * (i+1);
                continue;
            }

            matrix[i][j] = 0.1 * (i+1) * (j+1);
        }
    }

    return matrix;
}

static void WriteMatrix(double[][] matrix)
{
    foreach (var arr in matrix)
    {
        foreach (var item in arr)
        {
            Console.Write(double.Round(item, 2) + "\t");
        }

        Console.WriteLine();
    }
}