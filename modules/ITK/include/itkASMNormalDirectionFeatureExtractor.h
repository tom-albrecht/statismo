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

        typedef DiscreteGaussianImageFilter<TImage, TImage> GaussianFilterType;
        typedef BSplineInterpolateImageFunction<TImage, typename TImage::PixelType, typename TImage::PixelType> InterpolatedImageType;
        typedef typename InterpolatedImageType::Pointer InterpolatedImagePointerType;
        typedef typename InterpolatedImageType::CovariantVectorType InterpolatedGradientType;

    public:
        typedef typename statismo::Representer<TPointSet>::PointType PointType;
        typedef typename ASMNormalDirectionPointSampler<TPointSet>::Pointer SamplerPointerType;

        ASMNormalDirectionFeatureExtractor(SamplerPointerType sampler, float sigma): m_sampler(sampler), m_sigma(sigma) {
            std::cout << "FIXME: itk::ASMFE constructor" <<std::endl;
        }

        virtual ~ASMNormalDirectionFeatureExtractor() {
            std::cout << "FIXME: itk::ASMFE destructor" <<std::endl;
        }
        virtual statismo::VectorType ExtractFeatures(const TImage* const image, const TPointSet* const dataset, const PointType& point) const {
            statismo::VectorType features(m_sampler->GetNumberOfPoints());

            int i = 0;
            float sum = 0;

            InterpolatedImagePointerType interpolated = Interpolate(image);
            InterpolatedGradientType normal;

            normal.CastFrom(m_sampler->GetNormalForPoint(dataset, point));
            std::vector<PointType> samples = m_sampler->Sample(dataset, point);

            for (typename std::vector<PointType>::iterator it = samples.begin(); it != samples.end(); ++it) {
                PointType sample = *it;

                typename TImage::IndexType index;
                if (image->TransformPhysicalPointToIndex(sample, index)) {
                    features[i] = (float)(interpolated->EvaluateDerivative(sample) * normal);
                } else {
                    features[i] = 999;
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
            std::cout << "features for " <<point << " : " <<features <<std::endl;
            return features;
        }

    private:
        const SamplerPointerType m_sampler;
        const float m_sigma;


        InterpolatedImagePointerType Interpolate(const TImage* const image) const {
            InterpolatedImagePointerType inter = InterpolatedImageType::New();
            inter->SetSplineOrder(1);
            if (m_sigma != 0) {
                typename GaussianFilterType::Pointer smooth = GaussianFilterType::New();
                smooth->SetVariance(m_sigma * m_sigma);
                // FIXME: smooth->SetMaximumKernelWidth ???
                smooth->SetInput(image);
                smooth->Update();
                inter->SetInputImage(smooth->GetOutput());
            } else {
                inter->SetInputImage(image);
            }

            std::cout << "FIXME: blurring and interpolation done" << std::endl;
            return inter;
        }
    };

    template<typename TPointSet, typename TImage>
    class ASMNormalDirectionFeatureExtractorFactory : public statismo::ASMFeatureExtractorFactory<TPointSet, TImage> {

        typedef ASMNormalDirectionPointSampler<TPointSet> SamplerType;
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
