#include <stdio.h>
#include <math.h>

double var10(double x)
{
    return 2.0 * x * sin(x) - cos(x);
}

// в корне вычитаемое равно нулю, поэтому получается равенство x = phi10(x)
double phi10(double x)
{
    return x - var10(x) / 5.0;
}

double iteration(double (*phi)(double),
                 double a, double b,
                 double eps, double q,
                 int *steps, double *error)
{
    if (phi == NULL || steps == NULL || error == NULL)
    {
        return NAN;
    }

    *steps = 0;
    *error = NAN;

    if (!isfinite(a) || !isfinite(b) ||
        !isfinite(eps) || !isfinite(q) ||
        a > b || eps <= 0.0 || q < 0.0 || q >= 1.0)
    {
        return NAN;
    }

    double x = a / 2.0 + b / 2.0;

    // ограничиваем число подстановок на случай отсутствия сходимости
    while (*steps < 10000)
    {
        double next = phi(x);
        (*steps)++;

        // доказательство сходимости действует только на выбранном отрезке
        if (!isfinite(next) || next < a || next > b)
        {
            return NAN;
        }

        double delta = fabs(next - x);

        // по изменению приближений и коэффициенту сжатия
        // оцениваем расстояние от нового приближения до корня
        *error = q / (1.0 - q) * delta;

        // проверяем условие задания и оценку ошибки самого корня
        if (delta < eps && *error < eps)
        {
            return next;
        }

        x = next;
    }

    return NAN;
}

double bisection(double (*F)(double),
                double a, double b,
                double eps, int *steps)
{
    if (F == NULL || steps == NULL)
    {
        return NAN;
    }

    *steps = 0;

    if (!isfinite(a) || !isfinite(b) || !isfinite(eps) ||
        a > b || eps <= 0.0)
    {
        return NAN;
    }

    double fa = F(a);
    double fb = F(b);

    if (!isfinite(fa) || !isfinite(fb))
    {
        return NAN;
    }

    // корень на границе не требует дальнейшего сужения отрезка
    if (fa == 0.0)
    {
        return a;
    }

    if (fb == 0.0)
    {
        return b;
    }

    // без смены знака нельзя гарантировать наличие корня этим методом
    if ((fa > 0.0) == (fb > 0.0))
    {
        return NAN;
    }

    while (b - a >= eps)
    {
        double c = a / 2.0 + b / 2.0;

        // округление больше не позволяет сузить отрезок
        // хотя требуемое условие остановки ещё не выполнено
        if (c == a || c == b)
        {
            return NAN;
        }

        double fc = F(c);
        (*steps)++;

        if (!isfinite(fc))
        {
            return NAN;
        }

        if (fc == 0.0)
        {
            return c;
        }

        // сохраняем ту половину, на концах которой меняется знак
        if ((fa > 0.0) != (fc > 0.0))
        {
            b = c;
        }
        else
        {
            a = c;
            fa = fc;
        }
    }

    // ошибка середины не превышает половины длины итогового отрезка
    return a / 2.0 + b / 2.0;
}

int main(void)
{
    double a = 0.4;
    double b = 1.0;
    double eps = 1e-6;
    double table_root = 0.6533;

    int steps;
    double root = bisection(var10, a, b, eps, &steps);

    if (!isfinite(root))
    {
        fprintf(stderr, "Bisection failed!\n");
        return 1;
    }

    printf("Method: bisection\n");
    printf("Root: %.10f\n", root);
    printf("Table root: %.4f\n", table_root);
    printf("Steps: %d\n", steps);
    printf("Residual |F(x)|: %.3e\n", fabs(var10(root)));

    double iteration_error;
    int iteration_steps;

    double iteration_root = iteration(
        phi10, a, b, eps, 0.68,
        &iteration_steps, &iteration_error
    );

    if (!isfinite(iteration_root))
    {
        fprintf(stderr, "Iteration failed!\n");
        return 1;
    }

    printf("\nMethod: iteration\n");
    printf("Root: %.10f\n", iteration_root);
    printf("Table root: %.4f\n", table_root);
    printf("Error bound: %.3e\n", iteration_error);
    printf("Steps: %d\n", iteration_steps);
    printf("Residual |F(x)|: %.3e\n", fabs(var10(iteration_root)));

    return 0;
}