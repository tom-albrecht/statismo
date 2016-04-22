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


#ifndef STATISMO_MHFITTING_H
#define STATISMO_MHFITTING_H

#include "ActiveShapeModel.h"
#include "ASMFitting.h"
#include "ASMPointSampler.h"
#include <boost/thread.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/thread/future.hpp>

namespace statismo {



    class MHFittingConfiguration {
    private:

    public:
        MHFittingConfiguration(const ASMFittingConfiguration& asmFittingConfiguration) : m_asmFittingConfiguration(asmFittingConfiguration) { };
        const ASMFittingConfiguration& GetAsmFittingconfiguration() const {return m_asmFittingConfiguration; }

    private:
        ASMFittingConfiguration m_asmFittingConfiguration;
    };





    template<typename RigidTransformPointerType>
    class MHFittingResult {
    public:
        MHFittingResult() : m_isValid(false) {
        }

        MHFittingResult(bool isValid, void* chain, VectorType coefficients, RigidTransformPointerType rigidTransform) :
                m_isValid(isValid),
                m_chain(chain),
                m_coefficients(coefficients),
                m_rigidTransform(rigidTransform) { }

        bool IsValid() {
            return m_isValid;
        }

        void* GetChain() { return m_chain; }

        VectorType GetCoefficients() {
            return m_coefficients;
        }

        RigidTransformPointerType GetRigidTransform() {
            return m_rigidTransform;
        }

    private:
        bool m_isValid;
        void* m_chain;
        VectorType m_coefficients;
        RigidTransformPointerType m_rigidTransform;
    };




    template<typename TPointSet, typename TImage>
    class MHFittingStep {
        typedef ActiveShapeModel<TPointSet, TImage> ActiveShapeModelType;
        typedef ASMPointSampler<TPointSet, TImage> PointSamplerType;
        typedef typename Representer<TPointSet>::PointType PointType;
        typedef StatisticalModel<TPointSet> StatisticalModelType;
        typedef ASMFeatureExtractor<TPointSet, TImage> FeatureExtractorType;
        typedef ASMPreprocessedImage<TPointSet> PreprocessedImageType;
        typedef typename ActiveShapeModelType::RepresenterType::RigidTransformPointerType RigidTransformPointerType;
        typedef  ASMFittingStep<TPointSet, TImage> ASMFittingStepType;

        // FIXME change void*
        MHFittingStep(const MHFittingConfiguration &configuration, const ActiveShapeModelType *activeShapeModel,
                      const void* chain, const VectorType sourceCoefficients, const RigidTransformPointerType sourceTransform,
                       const PreprocessedImageType *targetImage, const std::vector<PointType> linePoints, const PointSamplerType *sampler)
                :
                m_configuration(configuration),
                m_model(activeShapeModel),
                m_sourceCoefficients(sourceCoefficients),
                m_sourceTransform(sourceTransform),
                m_target(targetImage),
                m_linePoints(linePoints),
                m_sampler(sampler) { }

    private:
        class ProfileResult {
        public:
            unsigned int pointId;
            PointType candidatePoint;
            PointType transformedCandidatePoint;
        };


    public:


        static MHFittingStep *Create(const MHFittingConfiguration &configuration,
                                      const ActiveShapeModelType *activeShapeModel,
                                      const void* chain,
                                      const VectorType sourceCoefficients,
                                      const RigidTransformPointerType sourceTransform,
                                      const PreprocessedImageType *targetImage,
                                      const std::vector<PointType> linePoints,
                                     const PointSamplerType *sampler) {
            return new MHFittingStep(configuration, activeShapeModel, chain, sourceCoefficients, sourceTransform, targetImage, linePoints, sampler);
        }

        ~MHFittingStep() {
        }

        MHFittingResult<RigidTransformPointerType> Perform() const {

            // runs a markov chain and returns only the accepted proposals

            ASMFittingStepType* asmFittingStep = ASMFittingStepType::Create(m_configuration.GetAsmFittingconfiguration(), m_model, m_sourceCoefficients, m_sourceTransform, m_target, m_sampler);
            ASMFittingResult<RigidTransformPointerType> result = asmFittingStep->Perform();
            MHFittingResult<RigidTransformPointerType> mhResult(true, m_chain, result.GetCoefficients(), result.GetRigidTransform());
            return mhResult;
        }

    private:
        const ActiveShapeModelType *m_model;
        void* m_chain;
        const VectorType m_sourceCoefficients;
        const RigidTransformPointerType m_sourceTransform;
        const PreprocessedImageType *m_target;
        std::vector<PointType> m_linePoints;
        const PointSamplerType *m_sampler;
        MHFittingConfiguration m_configuration;
    };
}

#endif //STATISMO_ASMFITTING_H

