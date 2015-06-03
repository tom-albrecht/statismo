//#include <itkImageFileReader.h>
//#include <itkActiveShapeModel.h>
//#include <itkStandardMeshRepresenter.h>
//#include <itkASMFitter.h>
//#include "itkASMTraits.h"
//#include "itkASMNormalDirectionGradientGaussianFeatureExtractor.h"
//typedef itk::ASMTraits<MeshType, ImageType> ASM;
//typedef itk::ImageFileReader<ASM::ImageType> ImageReaderType;
//typedef itk::ActiveShapeModel<ASM> ActiveShapeModelType;
//typedef itk::StandardMeshRepresenter<float, 3> RepresenterType;
//

#include <iostream>
#include <itkMesh.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkStandardMeshRepresenter.h>
#include <ASMFitter.h>
#include "itkASMNormalDirectionPointSampler.h"
#include "itkASMNormalDirectionFeatureExtractor.h"
#include "itkASMGaussianGradientImagePreprocessor.h"
#include "itkActiveShapeModel.h"
#include "itkASMFitter.h"

typedef itk::Mesh<float, 3> MeshType;
typedef itk::Image<float, 3> ImageType;
typedef itk::ActiveShapeModel<MeshType, ImageType> ActiveShapeModelType;
typedef itk::ImageFileReader<ImageType> ImageReaderType;

typedef itk::StandardMeshRepresenter<float, 3> RepresenterType;
typedef itk::StatisticalModel<MeshType> StatisticalModelType;
typedef vnl_vector<statismo::ScalarType> VnlVectorType;

VnlVectorType toVnlVector(const statismo::VectorType& v) {
    return VnlVectorType(v.data(), v.rows());

}

int main(int argc, char *argv[]) {

    std::cout << "Initializing..." << std::endl;
    // FIXME: these should go somewhere less "intrusive"
    statismo::ASMFeatureExtractorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMNormalDirectionFeatureExtractorFactory<MeshType, ImageType>::GetInstance());
    statismo::ASMImagePreprocessorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMGaussianGradientImagePreprocessorFactory<MeshType, ImageType>::GetInstance());

    std::string modelname("/tmp/asmbuild.h5");
    //std::string modelname("/tmp/asmLevel-1.h5");

    ActiveShapeModelType::Pointer aModel = ActiveShapeModelType::New();
    RepresenterType::Pointer representer = RepresenterType::New();
    aModel->Load(representer, modelname.c_str());

    itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::Pointer fitSampler = itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::New();
    fitSampler->SetNumberOfPoints(25);
    fitSampler->SetPointSpacing(0.5);

    statismo::ASMFitterConfiguration fitConfig(5,5,3);

    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName("/tmp/24.vtk");
    reader->Update();
    ImageType::Pointer image = reader->GetOutput();
    statismo::ASMPreprocessedImage<MeshType> *pimage = aModel->GetstatismoImplObj()->GetImagePreprocessor()->Preprocess(image);

    //FIXME: yuck. Essentially, I want a 0-initialized statismo vector of the correct size.
    //This should probably also be fixed to take/return itk Vectors, not statismo:: ones.
    StatisticalModelType::Pointer ssm = aModel->GetStatisticalModel();
    MeshType::Pointer mean = ssm->DrawMean();
    StatisticalModelType::VectorType mmean = ssm->ComputeCoefficientsForDataset(mean);
    statismo::VectorType smean = Eigen::Map<const statismo::VectorType>(mmean.data_block(), mmean.size());

    typedef itk::ASMFitterStep<MeshType, ImageType> FitterStepType;
    typedef typename FitterStepType::Pointer FitterStepPointerType;

    statismo::VectorType coeffs = smean;

    //FIXME: We'll need an actual "Fitter" type which does the iterations.
    FitterStepPointerType fitterStep = FitterStepType::New();
    fitterStep->SetModel(aModel);
    fitterStep->SetTarget(pimage);
    fitterStep->SetSampler(FitterStepType::SamplerPointerType(fitSampler.GetPointer()));
    fitterStep->SetConfiguration(fitConfig);

    std::cout << "Initialization done." << std::endl;

    for (int i =0; i < 25; ++i) {
        std::cout << "iteration: " << i << std::endl;
        fitterStep->SetCoefficients(coeffs);
        fitterStep->Update();
        statismo::ASMFitterResult result = fitterStep->GetOutput();
        if (!result.IsValid()) {
            std::cout << "invalid result, aborting " <<std::endl;
            exit(42);
        }
        coeffs = result.GetCoefficients();
        std::cout << "coeffs (adj)" << toVnlVector(coeffs) << std::endl;
    }
    return 0;
}
