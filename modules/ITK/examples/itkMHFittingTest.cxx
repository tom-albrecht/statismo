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
#include "MHFitting.h"
#include "itkMHFitting.h"
#include "itkMeshFileWriter.h"
#include "itkTimeProbe.h"
//#include "itkRigidTransformModelBuilder.h"


typedef itk::Mesh<float, 3> MeshType;
typedef itk::Image<float, 3> ImageType;
typedef itk::ActiveShapeModel<MeshType, ImageType> ActiveShapeModelType;
typedef itk::ImageFileReader<ImageType> ImageReaderType;

typedef itk::StandardMeshRepresenter<float, 3> RepresenterType;
typedef RepresenterType::RigidTransformType RigidTransformType;
typedef RepresenterType::PointType PointType;
typedef itk::StatisticalModel<MeshType> StatisticalModelType;
typedef itk::Euler3DTransform< float > TransformType;

typedef itk::MHFittingStep<MeshType, ImageType> FittingStepType;
typedef itk::MHFittingResult<MeshType, ImageType> FittingResultType;


// FIXME: these conversions have to go.
typedef vnl_vector<statismo::ScalarType> VnlVectorType;
VnlVectorType toVnlVector(const statismo::VectorType& v) {
    return VnlVectorType(v.data(), v.rows());
}
statismo::VectorType fromVnlVector(const VnlVectorType& v) {
    return Eigen::Map<const statismo::VectorType>(v.data_block(), v.size());

}

int main(int argc, char *argv[]) {


    std::cout << "Initializing..." << std::endl;
    // FIXME: these should go somewhere less "intrusive"
    statismo::ASMFeatureExtractorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMNormalDirectionFeatureExtractorFactory<MeshType, ImageType>::GetInstance());
    statismo::ASMImagePreprocessorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMGaussianGradientImagePreprocessorFactory<MeshType, ImageType>::GetInstance());

    std::string modelname("/export/skulls/data/shapes/ulna-right/aligned/registered-pami-ams/model-asm/asm-pca-3.h5");
//    std::string modelname("//home/marcel/data/ulna-right/test/asm-pca-3.h5");

    ActiveShapeModelType::Pointer aModel = ActiveShapeModelType::New();
    RepresenterType::Pointer representer = RepresenterType::New();
    aModel->Load(representer,  modelname.c_str());

    itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::Pointer fitSampler = itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::New();
    fitSampler->SetNumberOfPoints(25);
    fitSampler->SetPointSpacing(1);


    statismo::ASMFittingConfiguration asmFitConfig(3,5,3);
    statismo::MHFittingConfiguration mhFitConfig(asmFitConfig);



    // You should use here
    // > representer->ComputeRigidTransformFromLandmarks
    // if you have landmarks
    // Currently I only want an identity transform
    itk::VersorRigid3DTransform<float>::Pointer versorTransform = itk::VersorRigid3DTransform<float>::New();
    RigidTransformType::Pointer currentTransform(versorTransform.GetPointer());
    currentTransform->SetIdentity();

    // read and preprocess image
    ImageReaderType::Pointer reader = ImageReaderType::New();
    //reader->SetFileName("/export/skulls/data/shapes/submandibular_gland_l/aligned/initial/volume-ct/pddca-0522c0002.nii");
    reader->SetFileName("/export/skulls/data/shapes/ulna-right/aligned/initial/volume-ct/downsampled-2/vsd-0.nii");
    //reader->SetFileName("//home/marcel/data/ulna-right/test/image.nii");

    reader->Update();
    ImageType::Pointer image = reader->GetOutput();
    statismo::ASMPreprocessedImage<MeshType> *pimage = aModel->GetstatismoImplObj()->GetImagePreprocessor()->Preprocess(image);

    // just for testing
    aModel->SetStatisticalModel(aModel->GetStatisticalModel());
    statismo::VectorType coeffs = statismo::VectorType::Zero(aModel->GetStatisticalModel()->GetNumberOfPrincipalComponents());


    // a vector with all the point constraints that should be used within the fitting
    std::vector<PointType> linePoints;
    PointType t1,t2,t3,t4,t5,t6, t7, t8, t9;
    t1[0] = 63.4184455871582f; t1[1]=-56.53125f; t1[2] = 400.44671630859375f;
    t2[0] = 49.84854507446289f; t2[1]=-56.53125f; t2[2]=426.7882995605469f;
    t3[0] =-14.279170989990234f; t3[1]=-56.53125f; t3[2] =620.8013916015625f;
    t4[0]=-5.718493461608887; t4[1]=-56.53125f; t4[2]=635.908447265625f;
    t5[0]=10.8992919921875f; t5[1]=-56.53125f; t5[2]=623.8228149414062f;
    t6[0] = 49.68370056152344f; t6[1]=-51.84375f; t6[2]=510.5065002441406f;
    t7[0] =61.776329040527344f; t7[1]=-51.84375f; t7[2]=512.01806640625f;
    t8[0] = 35.07177734375f; t8[1]=-51.84375f, t8[2]=560.8924560546875f;
    t9[0] =46.15668487548828; t9[1]=-51.84375f; t9[2]=565.4271850585938f;


    linePoints.push_back(t1);
    linePoints.push_back(t2);
    linePoints.push_back(t3);
    linePoints.push_back(t4);
    linePoints.push_back(t5);
    linePoints.push_back(t6);
    linePoints.push_back(t7);
    linePoints.push_back(t8);
    linePoints.push_back(t9);




    // very ITK unlike, we use a init method instead of setting all fields manually.
    // This avoids 99% of all core dumps :-)
    FittingStepType::Pointer fittingStep = FittingStepType::New();
    fittingStep->init(pimage, linePoints, aModel, FittingStepType::SamplerPointerType(fitSampler.GetPointer()), mhFitConfig, currentTransform, coeffs);

    std::cout << "Initialization done." << std::endl;


    for (int i =1; i <= 1000; ++i) {

        std::cout << "iteration: " << i << std::endl;

        fittingStep->NextSample();
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

