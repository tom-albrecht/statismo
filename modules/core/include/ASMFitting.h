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


#ifndef STATISMO_ASMFITTING_H
#define STATISMO_ASMFITTING_H

#include "ActiveShapeModel.h"
#include "ASMPointSampler.h"

namespace statismo {

    class ASMFittingConfiguration {
    private:
        float m_featureDistanceThreshold;
        float m_pointDistanceThreshold;
        float m_modelCoefficientBounds;

    public:
        ASMFittingConfiguration(float featureDistanceThreshold, float pointDistanceThreshold, float modelCoefficientBounds) :
                m_featureDistanceThreshold(featureDistanceThreshold),
                m_pointDistanceThreshold(pointDistanceThreshold),
                m_modelCoefficientBounds(modelCoefficientBounds)
        {}


        float GetFeatureDistanceThreshold() const {
            return m_featureDistanceThreshold;
        }

        float GetPointDistanceThreshold() const {
            return m_pointDistanceThreshold;
        }

        float GetModelCoefficientBounds() const {
            return m_modelCoefficientBounds;
        }
    };

    template <typename RigidTransformPointerType>
    class ASMFittingResult {
    public:
        ASMFittingResult() : m_isValid(false) {
        }

        ASMFittingResult(bool isValid, VectorType coefficients, RigidTransformPointerType rigidTransform) :
            m_isValid(isValid),
            m_coefficients(coefficients),
            m_rigidTransform(rigidTransform) {}

        bool IsValid() {
            return m_isValid;
        }

         VectorType GetCoefficients() {
            return m_coefficients;
        }

        RigidTransformPointerType GetRigidTransform() {
            return m_rigidTransform;
        }

    private:
        bool m_isValid;
        VectorType m_coefficients;
        RigidTransformPointerType m_rigidTransform;
    };

    template <typename TPointSet, typename TImage>
    class ASMFittingStep {
        typedef ActiveShapeModel<TPointSet, TImage> ActiveShapeModelType;
        typedef ASMPointSampler<TPointSet, TImage> PointSamplerType;
        typedef typename Representer<TPointSet>::PointType PointType;
        typedef StatisticalModel<TPointSet> StatisticalModelType;
        typedef ASMFeatureExtractor<TPointSet, TImage> FeatureExtractorType;
        typedef ASMPreprocessedImage<TPointSet> PreprocessedImageType;
        typedef typename ActiveShapeModelType::RepresenterType::RigidTransformPointerType RigidTransformPointerType;

        ASMFittingStep(const ASMFittingConfiguration& configuration, const ActiveShapeModelType* activeShapeModel,
                  const VectorType sourceCoefficients, const RigidTransformPointerType sourceTransform, const PreprocessedImageType* targetImage, const PointSamplerType* sampler)
                :
                m_configuration(configuration),
                m_model(activeShapeModel),
                m_sourceCoefficients(sourceCoefficients),
                m_sourceTransform(sourceTransform),
                m_target(targetImage),
                m_sampler(sampler)
                //m_source(m_model->GetStatisticalModel()->DrawSample(sourceCoefficients, false))
        //m_featureExtractor(m_model->GetFeatureExtractor()->SetImage(m_targetImage)->SetPointset(m_inputMesh))
        {

        }

    public:

        static ASMFittingStep* Create(const ASMFittingConfiguration& configuration, const ActiveShapeModelType* activeShapeModel,
                                     const VectorType sourceCoefficients, const RigidTransformPointerType sourceTransform, const PreprocessedImageType* targetImage, const PointSamplerType* sampler) {
            return new ASMFittingStep(configuration, activeShapeModel, sourceCoefficients, sourceTransform, targetImage, sampler);
        }

        ~ASMFittingStep() {
            //FIXME: clean up m_source
        }

        ASMFittingResult<RigidTransformPointerType> Perform() const {
            typename StatisticalModelType::PointValueListType deplacements;
            std::vector<ASMProfile> profiles = m_model->GetProfiles();

            unsigned int dimensions = m_model->GetStatisticalModel()->GetRepresenter()->GetDimensions();

            std::vector<PointType> domainPoints = m_model->GetStatisticalModel()->GetRepresenter()->GetDomain().GetDomainPoints();

            PointSamplerType* sampler = m_sampler->CloneForTarget(m_model, m_sourceCoefficients, m_sourceTransform);
            FeatureExtractorType* fe = m_model->GetFeatureExtractor()->CloneForTarget(m_model, m_sourceCoefficients, m_sourceTransform);

            for (std::vector<ASMProfile>::const_iterator profile = profiles.begin(); profile != profiles.end(); ++profile) {
                unsigned pointId = (*profile).GetPointId();
                PointType profilePoint = m_model->GetStatisticalModel()->DrawSampleAtPoint(m_sourceCoefficients, pointId, false);
                PointType transformedProfilePoint = m_model->GetRepresenter()->TransformPoint(profilePoint, m_sourceTransform);
                PointType candidatePoint;
                float featureDistance = FindBestMatchingPointForProfile(candidatePoint, sampler, fe, transformedProfilePoint, (*profile).GetDistribution());
                if (featureDistance <= m_configuration.GetFeatureDistanceThreshold()) {
                    statismo::VectorType point(dimensions);
                    for (int i = 0; i < dimensions; ++i) {
                        point[i] = candidatePoint[i];
                    }
                    statismo::MultiVariateNormalDistribution marginal = m_model->GetMarginalAtPointId(pointId);
                    float pointDistance = marginal.MahalanobisDistance(point);

                    if (pointDistance <= m_configuration.GetPointDistanceThreshold()) {
                        PointType refPoint = domainPoints[pointId];
                        deplacements.push_back(std::make_pair(refPoint, candidatePoint));
                    } else {
                        // shape distance is larger than threshold
                    }
                } else {
                    // feature distance is larger than threshold
                }
            }

            sampler->Delete();
            fe->Delete();


            if (deplacements.size() > 0) {
                std::vector<PointType> reference;
                std::vector<PointType> target;

                for (typename StatisticalModelType::PointValueListType::const_iterator deplace = deplacements.begin(); deplace != deplacements.end(); ++deplace) {
                    reference.push_back(deplace->first);
                    target.push_back(deplace->second);
                }

                RigidTransformPointerType adjust = m_model->GetRepresenter()->ComputeRigidTransformFromLandmarks(reference, target);
                //RigidTransformPointerType adjust = m_model->GetRepresenter()->ComputeRigidTransformFromLandmarks(target, reference);
                if (adjust) {
                    std::cout << "constraints: " << deplacements.size() << std::endl;
//                    for(typename std::vector<PointType>::const_iterator pt = fixed.begin(); pt != fixed.end(); ++pt) {
//                        PointType ox = *pt;
//                        PointType tx = m_model->GetRepresenter()->TransformPoint(ox, adjust);
//                        std::cout << ox << " -> " << tx << std::endl;
//                    }
//                    std::cout << std::endl<< std::endl<< std::endl;
//                    for(typename std::vector<PointType>::const_iterator pt = moving.begin(); pt != moving.end(); ++pt) {
//                        PointType ox = *pt;
//                        PointType tx = m_model->GetRepresenter()->TransformPoint(ox, adjust);
//                        std::cout << ox << " -> " << tx << std::endl;
//                    }
//                    std::cout << std::endl<< std::endl<< std::endl;
                }
                if (adjust) {
                    // We need to transform back the candidate points to the model space to calculate coefficients
                    deplacements.clear();
                    int i = 0;
                    for(typename std::vector<PointType>::const_iterator ref = reference.begin(); ref != reference.end(); ++ref) {
                        PointType refPt = *ref;
                        PointType tgtPtT = target[i];
                        PointType tgtPtR = m_model->GetRepresenter()->TransformPoint(tgtPtT, adjust, true);
                        //std::cout << refPt << " -> " << tgtPtT << " -> " << tgtPtR <<std::endl;
                        deplacements.push_back(std::make_pair(refPt,tgtPtR));
                        ++i;
                    }

                }

                VectorType coeffs = m_model->GetStatisticalModel()->ComputeCoefficientsForPointValues(deplacements, 1e-6);

                float maxCoeff = m_configuration.GetModelCoefficientBounds();
                for (size_t i = 0; i < coeffs.size(); ++i) {
                    coeffs(i) = std::min(maxCoeff, std::max(-maxCoeff, coeffs(i)));
                }
                // , m_model->GetStatisticalModel()->DrawSample(coeffs)
                return ASMFittingResult<RigidTransformPointerType>(true, coeffs, adjust);

            } else {
                // Invalid result: coefficients with no content
                std::cout << "FIXME: returning invalid result" << std::endl;
                VectorType coeffs;
                return ASMFittingResult<RigidTransformPointerType>(false, coeffs, 0);
            }
        }
    private:
        const ActiveShapeModelType* m_model;
        const VectorType m_sourceCoefficients;
        const RigidTransformPointerType m_sourceTransform;
        const PreprocessedImageType* m_target;
        const PointSamplerType* m_sampler;
        ASMFittingConfiguration m_configuration;

        float FindBestMatchingPointForProfile(PointType &result, const PointSamplerType* sampler, const FeatureExtractorType* const fe, const PointType &profilePoint, const statismo::MultiVariateNormalDistribution &profile) const {
            float bestFeatureDistance = -1;

            std::vector<PointType> samples = sampler->Sample(profilePoint);

            for (typename std::vector<PointType>::const_iterator sample = samples.begin(); sample != samples.end(); ++sample) {
                statismo::VectorType features;
                bool ok = fe->ExtractFeatures(features, m_target, *sample);
                if (ok) {
                    float featureDistance = profile.MahalanobisDistance(features);
                    if (bestFeatureDistance == -1 || bestFeatureDistance > featureDistance) {
                        result = *sample;
                        bestFeatureDistance = featureDistance;
                    }
                }
            }
            return bestFeatureDistance;
        }
    };
}

#endif //STATISMO_ASMFITTING_H

