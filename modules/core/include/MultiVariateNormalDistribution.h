/*
 * This file is part of the statismo library.
 *
 * Author: Christoph Langguth (christoph.langguth@unibas.ch)
 *
 * Copyright (c) 2011-2015 University of Basel
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the project's author nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef STATISMO_MULTIVARIATENORMALDISTRIBUTION_H
#define STATISMO_MULTIVARIATENORMALDISTRIBUTION_H

#include "CommonTypes.h"
#include <cmath>
#include "StatismoUtils.h"

namespace statismo {
    class MultiVariateNormalDistribution {

    private:
        MatrixType covInv;

    public:

        VectorType mean;
        MatrixType covariance;


        MultiVariateNormalDistribution() {}

        MultiVariateNormalDistribution(VectorType _mean, MatrixType _covariance): mean(_mean), covariance(_covariance) {
            if (!Utils::PseudoInverse(covariance, covInv)) {
                //FIXME: throw some exception
            }
        }

        float logpdf(VectorType data) const {
            double mhDist = MahalanobisDistance(data);
            double detCov = 1;
            for (unsigned i = 0; i < covariance.rows(); ++i) {
                detCov *= covariance(i,i);
            }
            double normalizer = -0.5*mean.size() * std::log(2*M_PI) -0.5 * detCov;
            return -0.5 * mhDist + normalizer;
        }

        float MahalanobisDistance(VectorType data) const {
            VectorType x0 = data - mean;

            float d = sqrt(x0.dot((covInv * x0)));
            return d;
        }
    };
}
#endif //STATISMO_MULTIVARIATENORMALDISTRIBUTION_H
