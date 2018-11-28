/**
 * Copyright 2018 Brendan Duke.
 *
 * This file is part of Evaluating Derivatives.
 *
 * Evaluating Derivatives is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * Evaluating Derivatives is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Evaluating Derivatives. If not, see <http://www.gnu.org/licenses/>.
 */
#include <math.h>      // for pow
#include <stdint.h>    // for uint32_t
#include <stdio.h>     // for printf
#include <stdlib.h>    // for size_t, EXIT_SUCCESS

/**
 * array_size() - get the number of elements in array @arr.
 * @arr: array to be sized
 */
template<typename T, size_t N>
constexpr size_t
array_size(T (&)[N])
{
        return N;
}

static void
init_data(double *data, const size_t nelem, const double gamma)
{
        for (uint32_t i = 0;
             i < nelem;
             ++i) {
                data[i] = (i + 1)/gamma;
        }
}

static double
f(double *x, uint32_t n)
{
        double result = 0.0;

        for (uint32_t i = 0;
             i < n;
             ++i) {
                double x_i = x[i];
                result += x_i*x_i;
        }

        return result;
}

/**
 * Consider f(x) = \sum_{i = 1}^n x_i^2 with x_i = 1 for i = 1...n. Implement
 * in single and double precision.
 *
 * a) Examine errors [f(x + h*e_1) - f(x)]/h - 2 for n = 10^j and h = 10^{-k},
 * where 2 = 2x_1 is the first gradient component.
 *
 * Observe for which j and k the difference underflows to zero, and determine
 * the best possible approximation for j = 4.
 *
 * Check whether the order of summation or prescaling of the components by
 * \gamma, so that f(x) is calculated as \gamma^2 f(x/\gamma), makes any
 * difference.
 */
int main(void)
{
        constexpr size_t N = 1024;
        constexpr double gamma = 100.0;
        double x[N];

        init_data(x, array_size(x), gamma);

        double k = 0.0;
        for (;;) {
                double h = pow(10.0, -k);
                if (h == 0.0) {
                        printf("underflow for 10^-%.0f\n", k);
                        break;
                }

                for (uint32_t n = 1;
                     n < array_size(x);
                     n *= 10) {
                        x[0] = (1.0 + h)/gamma;
                        double err = f(x, n);
                        x[0] = 1.0/gamma;
                        err -= f(x, n);
                        err *= gamma*gamma;
                        if (err == 0.0) {
                                printf("difference underflown for k: %.0f n: %d\n", k, n);
                                if (n == 1)
                                        goto finish;
                                break;
                        }
                        err /= h;
                        err -= 2.0;

                        printf("k: %.0f n: %d err %f\n", k, n, err);
                }

                k += 1.0;
        }

finish:
        return EXIT_SUCCESS;
}
