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

/**
 * Consider f(x) = \sum_{i = 1}^n x_i^2 with x_i = 1 for i = 1...n. Implement
 * in single and double precision.
 *
 * a) Examine absolute errors [f(x + h*e_1) - f(x)]/h - 2 for n = 10^j and
 * h = 10^{-k}, where 2 = 2x_1 is the first gradient component.
 *
 * Observe for which j and k the difference underflows to zero, and determine
 * the best possible approximation for j = 4.
 *
 * Check whether the order of summation or prescaling of the components by
 * \gamma, so that f(x) is calculated as \gamma^2 f(x/\gamma), makes any
 * difference.
 */
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

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

/**
 * get_and_check_time_of_day() - Get time of day and store it in `tv`, checking
 * for errors.
 * @tv: Time of day output.
 */
static void
get_and_check_time_of_day(struct timeval *tv)
{
        int32_t status = gettimeofday(tv, NULL);
        assert(status == 0);
}

/**
 * get_seed_from_time_of_day() - Convenience function to return the 64-bit
 * micro-second part of the time of day.
 */
static uint64_t
get_seed_from_time_of_day(void)
{
        struct timeval seed;
        get_and_check_time_of_day(&seed);

        return seed.tv_usec;
}

/**
 * get_gsl_rng() - Allocates a GSL RNG and seeds it with the time in
 * microseconds.
 *
 * The caller owns the returned RNG.
 */
static gsl_rng *
get_gsl_rng(void)
{
        /**
         * NOTE(brendan): This must be done; passing the `gsl_rng_taus` pointer
         * directly to `gsl_rng_alloc` results in a SIGSEGV.
         */
        const gsl_rng_type *rng_type = gsl_rng_taus;
        gsl_rng *rng = gsl_rng_alloc(rng_type);
        assert(rng != NULL);

        gsl_rng_set(rng, get_seed_from_time_of_day());

        return rng;
}

/**
 * init_data_uniform() - Initializes `data` from a uniform distribution with
 * range [-a, a].
 * @data: The buffer to be initialized.
 * @rng: GSL RNG state to used to draw the random samples.
 * @nelem: Number of elements in data.
 * @a: Half-width of the uniform distribution.
 *
 * Returns the size in bytes of the entire `data` buffer initialized.
 */
template<typename T>
static void
init_data_uniform(T *data, gsl_rng *rng, const size_t nelem, const float a)
{
        for (uint32_t i = 0;
             i < nelem;
             ++i) {
                data[i] = gsl_ran_flat(rng, -a, a);
        }
}

template<typename T>
static T
f(T *x, uint32_t n)
{
        T result = 0.0;

        for (uint32_t i = 0;
             i < n;
             ++i) {
                T x_i = x[i];
                result += x_i*x_i;
        }

        return result;
}

int main(void)
{
        constexpr size_t N = 2048;
        constexpr uint32_t n = 12;
        float x_f[N];
        double x[N];

        gsl_rng *rng = get_gsl_rng();

        init_data_uniform(x_f, rng, array_size(x_f), 1.0f);
        init_data_uniform(x, rng, array_size(x), 1.0f);

        printf("%.5f\n", f(x_f, n));
        printf("%.5f\n", f(x, n));

        return EXIT_SUCCESS;
}
