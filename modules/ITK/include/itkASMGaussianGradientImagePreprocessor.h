#ifndef STATISMO_ITKASMGAUSSIANGRADIENTIMAGEPREPROCESSOR_H
#define STATISMO_ITKASMGAUSSIANGRADIENTIMAGEPREPROCESSOR_H

#include "ASMImagePreprocessor.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "HDF5Utils.h"

namespace itk {

    template <typename TPointSet, typename TImage>
    class ASMGaussianGradientPreprocessedImage: public statismo::ASMPreprocessedImage<TPointSet> {
    private:
        typedef DiscreteGaussianImageFilter<TImage, TImage> GaussianFilterType;
        typedef BSplineInterpolateImageFunction<TImage, typename TImage::PixelType, typename TImage::PixelType> InterpolatedImageType;
        typedef typename InterpolatedImageType::Pointer InterpolatedImagePointerType;
        //typedef typename InterpolatedImageType::CovariantVectorType CovariantVectorType;
        typedef CovariantVector<float, TImage::ImageDimension> CovariantVectorType;
        typedef vnl_vector<statismo::ScalarType> VectorType;


        const TImage* m_inputImage;
        const InterpolatedImagePointerType m_interpolatedImage;

        ASMGaussianGradientPreprocessedImage(const TImage* inputImage, const InterpolatedImagePointerType interpolatedImage): m_inputImage(inputImage), m_interpolatedImage(interpolatedImage) {}

        static statismo::VectorType fromVnlVector(const VectorType& v) {
            return Eigen::Map<const statismo::VectorType>(v.data_block(), v.size());

        }
    public:
        typedef typename statismo::Representer<TPointSet>::PointType PointType;

        virtual bool IsDefined(const PointType& point) const {
            typename TImage::IndexType index;
            // we don't care about the actual index, just whether it's present
            return m_inputImage->TransformPhysicalPointToIndex(point, index);
        }

        virtual statismo::VectorType Evaluate(const PointType& point) const {
            CovariantVectorType cv = m_interpolatedImage->EvaluateDerivative(point);
            return fromVnlVector(cv.GetVnlVector());
        }

        static ASMGaussianGradientPreprocessedImage* Create(const TImage* inputImage, const InterpolatedImagePointerType interpolatedImage) {
            return new ASMGaussianGradientPreprocessedImage(inputImage, interpolatedImage);
        }
    };

    template <typename TPointSet, typename TImage>
    class ASMGaussianGradientImagePreprocessor: public statismo::ASMImagePreprocessor<TPointSet, TImage> {
    private:
        typedef DiscreteGaussianImageFilter<TImage, TImage> GaussianFilterType;
        typedef BSplineInterpolateImageFunction<TImage, typename TImage::PixelType, typename TImage::PixelType> InterpolatedImageType;
        typedef typename InterpolatedImageType::Pointer InterpolatedImagePointerType;

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

            return inter;
        }

    public:
        typedef ASMGaussianGradientPreprocessedImage<TPointSet, TImage> PreprocessedImplType;

        ASMGaussianGradientImagePreprocessor(float sigma): m_sigma(sigma) {}

        virtual ASMGaussianGradientImagePreprocessor<TPointSet, TImage>* Clone() const {
            return new ASMGaussianGradientImagePreprocessor(m_sigma);
        };

        virtual PreprocessedImplType* Preprocess(const TImage* image) const {
            return PreprocessedImplType::Create(image, Interpolate(image));
        };
    };

    template<typename TPointSet, typename TImage>
    class ASMGaussianGradientImagePreprocessorFactory : public statismo::ASMImagePreprocessorFactory<TPointSet, TImage> {

        typedef ASMGaussianGradientImagePreprocessor<TPointSet, TImage> InstanceType;

    public:

        static const ASMGaussianGradientImagePreprocessorFactory *GetInstance() {
            static ASMGaussianGradientImagePreprocessorFactory *instance = new ASMGaussianGradientImagePreprocessorFactory();
            return instance;
        }

        virtual std::string GetDescriptor() const {
            return "builtin::GaussianGradient";
        }

        virtual const statismo::ASMImagePreprocessor<TPointSet, TImage> *Instantiate(
                const H5::Group &h5Group) const {

            float sigma = statismo::HDF5Utils::readFloat(h5Group, "stddev");
            InstanceType* instance = new InstanceType(sigma);
            return instance;
        }
    };


}
#endif //STATISMO_ITKASMGAUSSIANGRADIENTIMAGEPREPROCESSOR_H
