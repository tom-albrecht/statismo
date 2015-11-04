#include <iostream>
#include <itkMesh.h>
#include <itkImage.h>
#include <itkImageFileReader.h>

#include <itkStandardMeshRepresenter.h>
#include "ASMFitting.h"
#include <itkEuler3DTransform.h>
#include "itkASMNormalDirectionPointSampler.h"
#include "itkASMNormalDirectionFeatureExtractor.h"
#include "itkASMGaussianGradientImagePreprocessor.h"
#include "itkActiveShapeModel.h"
#include "itkASMFitting.h"
#include "itkMeshFileWriter.h"
//#include "itkRigidTransformModelBuilder.h"

typedef itk::Mesh<float, 3> MeshType;
typedef itk::Image<float, 3> ImageType;
typedef itk::ActiveShapeModel<MeshType, ImageType> ActiveShapeModelType;
typedef itk::ImageFileReader<ImageType> ImageReaderType;

typedef itk::StandardMeshRepresenter<float, 3> RepresenterType;
typedef RepresenterType::RigidTransformType RigidTransformType;
typedef itk::StatisticalModel<MeshType> StatisticalModelType;
typedef itk::Euler3DTransform< float > TransformType;

typedef itk::ASMFitting<MeshType, ImageType> FittingType;
typedef itk::ASMFittingStep<MeshType, ImageType> FittingStepType;
typedef itk::ASMFittingResult<MeshType, ImageType> FittingResultType;


// FIXME: these conversions have to go.
typedef vnl_vector<statismo::ScalarType> VnlVectorType;
VnlVectorType toVnlVector(const statismo::VectorType& v) {
    return VnlVectorType(v.data(), v.rows());
}
statismo::VectorType fromVnlVector(const VnlVectorType& v) {
    return Eigen::Map<const statismo::VectorType>(v.data_block(), v.size());

}


int main(int argc, char *argv[]) {


    /// TRANSFORM TESTING

    TransformType::Pointer transform = TransformType::New();
    TransformType::ParametersType parameters( 6 );
    // Euler angles (x,y,z)
    parameters[0] = 0.0; //0.1
    parameters[1] = 0.0; //0.2
    parameters[2] = 0.0; //0.3
    // Translation (x,y,z)
    parameters[3] = 0.0; //4.0
    parameters[4] = 0.0; //5.0
    parameters[5] = 0.0; //6.0
    transform->SetParameters( parameters );
#if defined(ITK_FIXED_PARAMETERS_ARE_DOUBLE) // After 4.8.1
    TransformType::FixedParametersType fixedParameters(3);
#else                                         //Pre 4.8.1
    TransformType::ParametersType fixedParameters(3);
#endif
    // rotation center
    fixedParameters[0] = -3.5;
    fixedParameters[1] = -4.5;
    fixedParameters[2] = -5.5;
    transform->SetFixedParameters( fixedParameters );
    std::cout << "Original transform: " << transform << std::endl;

    // TRANSFORM

    if (false) {
        std::cout << "Exiting." << std::endl;
        return 0;
    }


    std::cout << "Initializing..." << std::endl;
    // FIXME: these should go somewhere less "intrusive"
    statismo::ASMFeatureExtractorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMNormalDirectionFeatureExtractorFactory<MeshType, ImageType>::GetInstance());
    statismo::ASMImagePreprocessorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMGaussianGradientImagePreprocessorFactory<MeshType, ImageType>::GetInstance());

    std::string modelname("/export/skulls/data/shapes/optic_nerve_l/aligned/registered-0002/model-asm/asm-5.h5");
    ActiveShapeModelType::Pointer aModel = ActiveShapeModelType::New();
    RepresenterType::Pointer representer = RepresenterType::New();
    aModel->Load(representer,  modelname.c_str());

    itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::Pointer fitSampler = itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::New();
    fitSampler->SetNumberOfPoints(25);
    fitSampler->SetPointSpacing(1);

    statismo::ASMFittingConfiguration fitConfig(5,5,3);

    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName("/tmp/onl.vtk");
    reader->Update();
    ImageType::Pointer image = reader->GetOutput();
    statismo::ASMPreprocessedImage<MeshType> *pimage = aModel->GetstatismoImplObj()->GetImagePreprocessor()->Preprocess(image);

    statismo::VectorType coeffs = statismo::VectorType::Zero(aModel->GetStatisticalModel()->GetNumberOfPrincipalComponents());
    RigidTransformType::Pointer currentTransform = RigidTransformType::Pointer(transform.GetPointer());



//    //FIXME: To be implemented
//    FittingType::Pointer fitting = FittingType::New();
//    fitting->SetModel(aModel);
//    fitting->SetTarget(pimage);
//    fitting->SetSampler(FittingStepType::SamplerPointerType(fitSampler.GetPointer()));
//    fitting->SetConfiguration(fitConfig);
//    fitting->SetInitialRigidTransformation(RigidTransformType::Pointer(transform.GetPointer()));
//    fitting->SetInitialCoefficients(coeffs);
//    fitting->SetNumberOfIterations(25);
//    fitting->Update();



    FittingStepType::Pointer fittingStep = FittingStepType::New();
    fittingStep->SetModel(aModel);
    fittingStep->SetTarget(pimage);
    fittingStep->SetSampler(FittingStepType::SamplerPointerType(fitSampler.GetPointer()));
    fittingStep->SetConfiguration(fitConfig);

    std::cout << "Initialization done." << std::endl;

    for (int i =0; i < 25; ++i) {
        std::cout << "iteration: " << i << std::endl;
        fittingStep->SetCoefficients(coeffs);
        fittingStep->SetRigidTransformation(currentTransform);
        fittingStep->Update();
        FittingResultType::Pointer result = fittingStep->GetOutput();
        if (!result->IsValid()) {
            std::cout << "invalid result, aborting " <<std::endl;
            exit(42);
        }
        coeffs = fromVnlVector(result->GetCoefficients());
        currentTransform = result->GetRigidTransformation();
        std::cout << "coeffs (adj)" << toVnlVector(coeffs) << std::endl;
        if (currentTransform) {
            std::cout << "Writing result of iteration " << i << std::endl;
            itk::MeshFileWriter<MeshType>::Pointer writer = itk::MeshFileWriter<MeshType>::New();
            std::stringstream filename;
            filename << "/tmp/itkmesh-" << i << ".vtk";
            writer->SetFileName(filename.str());
            MeshType::Pointer ms = result->GetMesh();
            writer->SetInput(ms);
            writer->Update();
        }
    }
    return 0;
}
