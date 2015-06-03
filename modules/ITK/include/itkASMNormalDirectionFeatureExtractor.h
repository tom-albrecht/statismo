///*
// * This file is part of the statismo library.
// *
// * Author: Christoph Langguth (christoph.langguth@unibas.ch)
// *
// * Copyright (c) 2011-2015 University of Basel
// * All rights reserved.
// *
// * Redistribution and use in source and binary forms, with or without
// * modification, are permitted provided that the following conditions
// * are met:
// *
// * Redistributions of source code must retain the above copyright notice,
// * this list of conditions and the following disclaimer.
// *
// * Redistributions in binary form must reproduce the above copyright
// * notice, this list of conditions and the following disclaimer in the
// * documentation and/or other materials provided with the distribution.
// *
// * Neither the name of the project's author nor the names of its
// * contributors may be used to endorse or promote products derived from
// * this software without specific prior written permission.
// *
// * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// *
// */
//
#ifndef STATISMO_itkAsmNormalDirectionFeatureExtractor_H
#define STATISMO_itkAsmNormalDirectionFeatureExtractor_H

#include "ASMFeatureExtractor.h"
#include "itkASMFeatureExtractor.h"
#include "HDF5Utils.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkASMNormalDirectionPointSampler.h"
#include "itkDiscreteGaussianImageFilter.h"

namespace itk {
    template<typename TPointSet, typename TImage>
    class ASMNormalDirectionFeatureExtractor : public statismo::ASMFeatureExtractor<TPointSet, TImage> {
    private:
        typedef statismo::ASMPreprocessedImage<TPointSet> PreprocessedImageType;
        typedef statismo::ActiveShapeModel<TPointSet, TImage> ActiveShapeModelType;
        typedef statismo::VectorType VectorType;
        typedef CovariantVector<float, TImage::ImageDimension> CovariantVectorType;
        typedef vnl_vector<statismo::ScalarType> VnlVectorType;
        typedef typename ASMNormalDirectionPointSampler<TPointSet, TImage>::Pointer SamplerPointerType;

        const SamplerPointerType m_sampler;
        const float m_sigma;

        static statismo::VectorType fromVnlVector(const VnlVectorType& v) {
            return Eigen::Map<const statismo::VectorType>(v.data_block(), v.size());
        }

    public:
        typedef typename statismo::Representer<TPointSet>::PointType PointType;

        ASMNormalDirectionFeatureExtractor(SamplerPointerType sampler, float sigma): m_sampler(sampler), m_sigma(sigma) {
        }

        virtual ~ASMNormalDirectionFeatureExtractor() {
        }

        virtual void Delete() {
            m_sampler->Delete();
            delete this;
        }

        virtual ASMNormalDirectionFeatureExtractor* CloneForTarget(const ActiveShapeModelType* const model, const VectorType& coefficients) const {
            SamplerPointerType copySampler = m_sampler->CloneForTarget(model, coefficients);
            return new ASMNormalDirectionFeatureExtractor(copySampler, m_sigma);
        }


        virtual bool ExtractFeatures(statismo::VectorType& output, const PreprocessedImageType* const image, const PointType& point) const {

            statismo::VectorType features(m_sampler->GetNumberOfPoints());

            int i = 0;
            float sum = 0;

            CovariantVectorType normal;
            normal.CastFrom(m_sampler->GetNormalForPoint(point));

            std::vector<PointType> samples = m_sampler->Sample(point);
            for (typename std::vector<PointType>::iterator it = samples.begin(); it != samples.end(); ++it) {
                PointType sample = *it;

                typename TImage::IndexType index;
                if (image->IsDefined(sample)) {
                    //features[i] = (float)(interpolated->EvaluateDerivative(sample) * normal);
                    VectorType gradient = image->Evaluate(sample);
                    VectorType snormal = fromVnlVector(normal.GetVnlVector());
                    float directionGradient = gradient.dot(snormal);
                    features[i] = directionGradient;
                } else {
                    // fail fast
                    //std::cout << "point " << point << " is out of bounds, returning empty feature vector" << std::endl;
                    return false;
                }
                // keep track of sum of absolute values
                if (features[i] < 0) {
                    sum -= features[i];
                } else {
                    sum += features[i];
                }
                ++i;
            }
            //normalize
            if (sum > 0) {
                for (i = 0; i < features.size(); ++i) {
                    features[i] /= sum;
                }
            }
            output = features;
            return true;
        }
    };

    template<typename TPointSet, typename TImage>
    class ASMNormalDirectionFeatureExtractorFactory : public statismo::ASMFeatureExtractorFactory<TPointSet, TImage> {

        typedef ASMNormalDirectionPointSampler<TPointSet, TImage> SamplerType;
        typedef typename SamplerType::Pointer SamplerPointerType;
        typedef ASMNormalDirectionFeatureExtractor<TPointSet, TImage> InstanceType;

    public:

        static const ASMNormalDirectionFeatureExtractorFactory *GetInstance() {
            static ASMNormalDirectionFeatureExtractorFactory *instance = new ASMNormalDirectionFeatureExtractorFactory();
            return instance;
        }

        virtual std::string GetDescriptor() const {
            return "builtin::NormalDirection";
        }

        virtual const statismo::ASMFeatureExtractor<TPointSet, TImage> *Instantiate(
                const H5::Group &h5Group) const {

            unsigned int numberOfPoints = statismo::HDF5Utils::readInt(h5Group, "numberOfPoints");
            float spacing = statismo::HDF5Utils::readFloat(h5Group, "spacing");
            SamplerPointerType sampler = SamplerType::New();
            sampler->SetNumberOfPoints(numberOfPoints);
            sampler->SetPointSpacing(spacing);

            InstanceType* instance = new InstanceType(sampler, 0);
            return instance;
            //m_sigma = statismo::HDF5Utils::readFloat(h5Group, "sigma");;
        }
    };

}
#endif //STATISMO_itkAsmNormalDirectionFeatureExtractor_H
