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


#ifndef STATISMO_ASMFITTER_H
#define STATISMO_ASMFITTER_H

#include "ActiveShapeModel.h"
#include "ASMPointSampler.h"

namespace statismo {

    class ASMFitterConfiguration {
    private:
        float m_maxFeatureDistance;
        float m_maxShapeDistance;
        float m_maxCoefficientDistance;

    public:
        ASMFitterConfiguration(float maxFeatureDistance, float maxShapeDistance, float maxCoefficientDistance) :
                m_maxFeatureDistance(maxFeatureDistance),
                m_maxShapeDistance(maxShapeDistance),
                m_maxCoefficientDistance(maxCoefficientDistance)
        {}


        float GetMaxFeatureDistance() const {
            return m_maxFeatureDistance;
        }

        float GetMaxShapeDistance() const {
            return m_maxShapeDistance;
        }

        float GetMaxCoefficientDistance() const {
            return m_maxCoefficientDistance;
        }
    };

    template <typename TPointSet, typename TImage>
    class ASMFitterResult {
    public:
        ASMFitterResult(bool isValid, VectorType coefficients, TPointSet* mesh) :
            m_isValid(isValid),
            m_coefficients(coefficients),
            m_mesh(mesh)
                    {
        }

        bool IsValid() {
            return m_isValid;
        }

         VectorType GetCoefficients() {
            return m_coefficients;
        }

         TPointSet* GetMesh() {
            return m_mesh;
        }

    private:
        bool m_isValid;
        VectorType m_coefficients;
        TPointSet* m_mesh;
    };

    template <typename TPointSet, typename TImage>
    class ASMFitter {
        typedef ASMFitter<TPointSet, TImage> ASMFitterType;
        typedef ActiveShapeModel<TPointSet, TImage> ActiveShapeModelType;
        typedef ASMPointSampler<TPointSet> PointSamplerType;
        typedef typename Representer<TPointSet>::PointType PointType;
        typedef StatisticalModel<TPointSet> StatisticalModelType;
        typedef ASMFeatureExtractor<TPointSet, TImage> FeatureExtractorType;


        ASMFitter(const ASMFitterConfiguration& configuration, const ActiveShapeModelType* activeShapeModel,
                  const TPointSet* inputMesh, const TImage* targetImage, const PointSamplerType* sampler)
                :
                m_configuration(configuration),
                m_model(activeShapeModel),
                m_inputMesh(inputMesh),
                m_targetImage(targetImage),
                m_sampler(sampler),
                m_featureExtractor(m_model->GetFeatureExtractor()->SetImage(m_targetImage)->SetPointset(m_inputMesh))
        {

        }

    public:

        static ASMFitter<TPointSet, TImage>* Create(const ASMFitterConfiguration& configuration, const ActiveShapeModelType* activeShapeModel,
                                                    const TPointSet* inputMesh, const TImage* targetImage, const PointSamplerType* sampler) {
            return new ASMFitter(configuration, activeShapeModel, inputMesh, targetImage, sampler);
        }


        ASMFitterResult<TPointSet, TImage> Fit() {
            typename StatisticalModelType::PointValueListType constraints;
            std::vector<ASMProfile> profiles = m_model->GetProfiles();

            unsigned int dimensions = m_model->GetStatisticalModel()->GetRepresenter()->GetDimensions();

            std::vector<PointType> domainPoints = m_model->GetStatisticalModel()->GetRepresenter()->GetDomain().GetDomainPoints();

            for (std::vector<ASMProfile>::const_iterator profile = profiles.begin(); profile != profiles.end(); ++profile) {
                unsigned pointId = (*profile).GetPointId();
                //std::cout << "evaluating @ pointId " << pointId << std::endl;
                PointType transformedPoint;
                float featureDistance = FindBestMatchingPointForProfile(transformedPoint, pointId, (*profile).GetDistribution());
                if (m_configuration.GetMaxFeatureDistance() == 0 || featureDistance <= m_configuration.GetMaxFeatureDistance()) {
                    statismo::VectorType point(dimensions);
                    for (int i = 0; i < dimensions; ++i) {
                        point[i] = transformedPoint[i];
                    }
                    statismo::MultiVariateNormalDistribution marginal = m_model->GetMarginalAtPointId(pointId);
                    float pointDistance = marginal.MahalanobisDistance(point);

                    if (m_configuration.GetMaxShapeDistance() == 0 || pointDistance <= m_configuration.GetMaxShapeDistance()) {
                        PointType refPoint = domainPoints[pointId];
                        constraints.push_back(std::make_pair(refPoint, transformedPoint));
                    } else {
                        // shape distance is larger than threshold
                    }
                } else {
                    // feature distance is larger than threshold
                }
            }


            if (constraints.size() > 0) {
                //FIXME: find rigid transformation which minimizes the distances between the constraint pair points

                VectorType coeffs = m_model->GetStatisticalModel()->ComputeCoefficientsForPointValues(constraints, 1e-6);

                float maxCoeff = m_configuration.GetMaxCoefficientDistance();
                if (maxCoeff > 0) {
                    for (size_t i = 0; i < coeffs.size(); ++i) {
                        coeffs(i) = std::min(maxCoeff, std::max(-maxCoeff, coeffs(i)));
                    }
                }
                return ASMFitterResult<TPointSet, TImage>(true, coeffs, m_model->GetStatisticalModel()->DrawSample(coeffs));

            } else {
                // Invalid result: coefficients with no content, mesh is null
                VectorType coeffs;
                return ASMFitterResult<TPointSet, TImage>(false, coeffs , 0);
            }

        }

        virtual ASMFitterType* SetMesh(TPointSet* mesh) {
               return  new ASMFitter(m_configuration, m_model, mesh, m_targetImage, m_sampler);
        };

    protected:
        ASMFitterConfiguration m_configuration;
        const ActiveShapeModelType* m_model;
        const TPointSet* m_inputMesh;
        const TImage* m_targetImage;
        const PointSamplerType* m_sampler;
        const FeatureExtractorType* m_featureExtractor;


        //virtual typename ASM::FitterResultPointerType InstantiateResult() = 0;


    private:
        float FindBestMatchingPointForProfile(PointType &result, unsigned pointId, const statismo::MultiVariateNormalDistribution &profile) {
            float bestFeatureDistance = -1;

            std::vector<PointType> samples = m_sampler->SampleAtPointId(pointId);

            for (typename std::vector<PointType>::const_iterator sample = samples.begin(); sample != samples.end(); ++sample) {
                statismo::VectorType features = m_featureExtractor->ExtractFeatures(*sample);
                float featureDistance = profile.MahalanobisDistance(features);
                if (bestFeatureDistance == -1 || bestFeatureDistance > featureDistance) {
                    result = *sample;
                    bestFeatureDistance = featureDistance;
                }
            }
            return bestFeatureDistance;
        }
    };
}

#endif //STATISMO_ASMFITTER_H

