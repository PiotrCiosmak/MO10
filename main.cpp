#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

constexpr int N{1000};

double analyticalFormula(double t)
{
    return 1 - exp(-10.0 * (t + atan(t)));
}

double directEulerMethod(double step)
{
    double exactValue, error{}, t{step}, y{};//warunek poczatkowy y(0) = 0

    for (int i = 0; i < N; ++i)
    {
        y = y + step * (-((10.0 * t * t + 20.0) / (t * t + 1.0)) * (y - 1.0));
        exactValue = fabs(analyticalFormula(t) - y);
        if (exactValue > error)
            error = exactValue;
        t += step;
    }
    return error;
}

double indirectEulerMethod(double step)
{
    double exactValue, fraction, error = 0.0, t = step, y = 0.0;//warunek poczatkowy y(0) = 0

    for (int i = 0; i < N; ++i)
    {
        fraction = (10.0 * (t + step) * (t + step) + 20.0) / ((t + step) * (t + step) + 1.0);
        y = (y + step * fraction) / (1 + step * fraction);
        exactValue = fabs(analyticalFormula(t) - y);
        if (exactValue > error)
            error = exactValue;
        t += step;
    }
    return error;
}

double trapezoidalMethod(double step)
{
    double exactValue, fractionN, fractionNPlus, error = 0.0, t = step, y = 0.0;//warunek poczatkowy y(0) = 0

    for (int i = 0; i < N; ++i)
    {
        fractionN = ((10.0 * t * t + 20.0) / (t * t + 1.0));
        fractionNPlus = (10.0 * (t + step) * (t + step) + 20.0) / ((t + step) * (t + step) + 1.0);;
        y = ((-step / 2.0) * (fractionN * (y - 1.0) - fractionNPlus) + y) / (1.0 + (step / 2.0) * fractionNPlus);
        exactValue = fabs(analyticalFormula(t) - y);
        if (exactValue > error)
            error = exactValue;
        t += step;
    }
    return error;

}

double directEulerMethodGraph(double step, double tmax)
{
    double y = 0.0; //warunek poczatkowy y(0) = 0

    for (double i = 0.0; i < tmax; i += step)
        y = y + step * (-((10.0 * i * i + 20.0) / (i * i + 1.0)) * (y - 1.0));
    return y;
}

double indirectEulerMethodGraph(double step, double tmax)
{
    double fraction, y = 0.0; //warunek poczatkowy y(0) = 0

    for (double i = 0.0; i < tmax; i += step)
    {
        fraction = (10.0 * (i + step) * (i + step) + 20.0) / ((i + step) * (i + step) + 1.0);
        y = (y + step * fraction) / (1 + step * fraction);
    }
    return y;
}

double trapezoidalMethodGraph(double step, double tmax)
{
    double fractionN, fractionNPlus, y = 0.0; //warunek poczatkowy y(0) = 0

    for (double i = 0.0; i < tmax; i += step)
    {
        fractionN = ((10.0 * i * i + 20.0) / (i * i + 1.0));
        fractionNPlus = (10.0 * (i + step) * (i + step) + 20.0) / ((i + step) * (i + step) + 1.0);;
        y = ((-step / 2.0) * (fractionN * (y - 1.0) - fractionNPlus) + y) / (1.0 + (step / 2.0) * fractionNPlus);
    }
    return y;
}

int main()
{
    fstream directEulerMethodErrors("directEulerMethodErrors.txt", fstream::out);
    fstream indirectEulerMethodErrors("indirectEulerMethodErrors.txt", fstream::out);
    fstream trapezoidalMethodErrors("trapezoidalMethodErrors.txt", fstream::out);
    for (double step = 0.1; step > 1e-20; step /= 2)
    {
        directEulerMethodErrors << log10(step) << "," << log10(directEulerMethod(step)) << endl;
        indirectEulerMethodErrors << log10(step) << "," << log10(indirectEulerMethod(step)) << endl;
        trapezoidalMethodErrors << log10(step) << "," << log10(trapezoidalMethod(step)) << endl;
    }
    directEulerMethodErrors.close();
    indirectEulerMethodErrors.close();
    trapezoidalMethodErrors.close();


    fstream directEulerMethodStable("directEulerMethodStable.txt", fstream::out);
    fstream indirectEulerMethod("indirectEulerMethod.txt", fstream::out);
    fstream trapezoidalMethod("trapezoidalMethod.txt", fstream::out);
    for (double t = 0; t < 5; t += 0.01)
    {
        anl = analyticalFormula(t);
        mpe = indirectEulerMethodGraph(0.01, t);
        mt = trapezoidalMethodGraph(0.01, t);
        mbe_s = directEulerMethodGraph(0.01, t);
        //dane_1 << t << " " << anl << " " << mpe << " " << mt << " " << mbe_s << endl;
    }
    directEulerMethodStable.close();
    indirectEulerMethod.close();
    trapezoidalMethod.close();


    fstream directEulerMethodUnStable("directEulerMethodUnStable.txt", fstream::out);
    for (double t = 0; t < 5; t += 0.15)
        directEulerMethodUnStable << t << " " << directEulerMethodGraph(0.15, t) << endl;
    directEulerMethodUnStable.close();
}