#include <stdio.h>
#include <math.h>

double var10(double x)
{
    return 2.0 * x * sin(x) - cos(x);
}

double derivative10(double x)
{
    return 3.0 * sin(x) + 2.0 * x * cos(x);
}

// в корне вычитаемое равно нулю, поэтому получается равенство x = phi10(x)
double phi10(double x)
{
    return x - var10(x) / 5.0;
}

double iteration(double (*phi)(double),
                 double a, double b,
                 double eps, int *steps)
{
    if (phi == NULL || steps == NULL)
    {
        return NAN;
    }

    *steps = 0;

    if (!isfinite(a) || !isfinite(b) ||
        !isfinite(eps) || a > b || eps <= 0.0)
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

        // Останавливаемся, когда соседние приближения отличаются меньше eps.
        if (delta < eps)
        {
            return next;
        }

        x = next;
    }

    return NAN;
}

double dichotomition(double (*F)(double),
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

double newton(double (*F)(double), double (*dF)(double),
              double a, double b,
              double eps, int *steps)
{
    if (F == NULL || dF == NULL || steps == NULL)
    {
        return NAN;
    }

    *steps = 0;

    if (!isfinite(a) || !isfinite(b) || !isfinite(eps) ||
        a > b || eps <= 0.0)
    {
        return NAN;
    }

    double x = a / 2.0 + b / 2.0;

    while (*steps < 10000)
    {
        double fx = F(x);
        if (!isfinite(fx))
        {
            return NAN;
        }

        // Если уже попали в корень, поправка больше не нужна.
        if (fx == 0.0)
        {
            return x;
        }

        double dfx = dF(x);
        if (!isfinite(dfx) || dfx == 0.0)
        {
            return NAN;
        }

        // Пересечение касательной с осью x даёт новое приближение.
        double next = x - fx / dfx;
        (*steps)++;

        // При неудачном шаге прекращаем поиск на выбранном отрезке.
        if (!isfinite(next) || next < a || next > b)
        {
            return NAN;
        }

        if (fabs(next - x) < eps)
        {
            return next;
        }

        x = next;
    }

    return NAN;
}

int main(void)
{
    double a = 0.4;
    double b = 1.0;
    double eps = 1e-6;
    double table_root = 0.6533;

    int steps;
    double root = dichotomition(var10, a, b, eps, &steps);

    if (!isfinite(root))
    {
        fprintf(stderr, "dichotomition failed!\n");
        return 1;
    }

    printf("Method: dichotomition\n");
    printf("Root: %.10f\n", root);
    printf("Table root: %.4f\n", table_root);
    printf("Steps: %d\n", steps);
    printf("Residual |F(x)|: %.3e\n", fabs(var10(root)));

    int iteration_steps;

    double iteration_root = iteration(
        phi10, a, b, eps, &iteration_steps
    );

    if (!isfinite(iteration_root))
    {
        fprintf(stderr, "Iteration failed!\n");
        return 1;
    }

    printf("\nMethod: iteration\n");
    printf("Root: %.10f\n", iteration_root);
    printf("Table root: %.4f\n", table_root);
    printf("Steps: %d\n", iteration_steps);
    printf("Residual |F(x)|: %.3e\n", fabs(var10(iteration_root)));

    int newton_steps;
    double newton_root = newton(var10, derivative10, a, b, eps, &newton_steps);

    if (!isfinite(newton_root))
    {
        return 1;
    }

    printf("\nMethod: Newton\n");
    printf("Root: %.10f\n", newton_root);
    printf("Table root: %.4f\n", table_root);
    printf("Steps: %d\n", newton_steps);
    printf("Residual |F(x)|: %.3e\n", fabs(var10(newton_root)));

    return 0;
}
