#include <iostream>
#include <itkMesh.h>
#include <itkImage.h>
#include <itkImageFileReader.h>

#include <itkStandardMeshRepresenter.h>
#include "ASMFitting.h"
#include <itkEuler3DTransform.h>
#include <ITK-4.5/itkCenteredTransformInitializer.h>
#include <ITK-4.5/itkCenteredVersorTransformInitializer.h>
#include "itkASMNormalDirectionPointSampler.h"
#include "itkASMNormalDirectionFeatureExtractor.h"
#include "itkASMGaussianGradientImagePreprocessor.h"
#include "itkActiveShapeModel.h"
#include "itkASMFitting.h"
#include "MHFitting.h"
#include "itkMHFitting.h"
#include "itkMeshFileWriter.h"
#include "itkTimeProbe.h"
#include "itkMeshFileReader.h"
#include "../cli/utils/statismo-fitting-utils.h"
#include "itkReducedVarianceModelBuilder.h"
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

typedef itk::MHFittingStepper<MeshType, ImageType> FittingStepType;
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

    //std::string modelname("/export/skulls/data/shapes/ulna-right/aligned/registered-pami-ams/model-asm/asm-pca-3.h5");
    std::string modelname("//tmp/fancyasm.h5");
    //std::string modelname("//home/marcel/data/ulna-right/test/asm-pca-3.h5");

    ActiveShapeModelType::Pointer aModel = ActiveShapeModelType::New();
    RepresenterType::Pointer representer = RepresenterType::New();
    aModel->Load(representer,  modelname.c_str());

    itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::Pointer fitSampler = itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::New();
    fitSampler->SetNumberOfPoints(25);
    fitSampler->SetPointSpacing(1);


    statismo::ASMFittingConfiguration asmFitConfig(3,5,3);
    statismo::MHFittingConfiguration mhFitConfig(asmFitConfig);

    // read and preprocess image
    ImageReaderType::Pointer reader = ImageReaderType::New();
    //reader->SetFileName("/export/skulls/data/shapes/submandibular_gland_l/aligned/initial/volume-ct/pddca-0522c0002.nii");
    //reader->SetFileName("/export/skulls/data/shapes/ulna-right/aligned/initial/volume-ct/downsampled-2/vsd-0.nii");
    reader->SetFileName("/export/skulls/data/shapes/esophagus/aligned/initial/volume-ct/varian-0021.nii");
//    reader->SetFileName("//home/marcel/data/ulna-right/test/image.nii");

    reader->Update();
    ImageType::Pointer image = reader->GetOutput();


    // You should use here
    // > representer->ComputeRigidTransformFromLandmarks
    // if you have landmarks
    // Currently I only want an identity transform
    itk::VersorRigid3DTransform<float>::Pointer currentTransform = itk::VersorRigid3DTransform<float>::New();
    //RigidTransformType::Pointer currentTransform(versorTransform.GetPointer());
    currentTransform->SetIdentity();
    itk::Point<float, 3> center;
    for (unsigned d =0; d < 3; ++d) {
        center[d] = image->GetOrigin()[d] + image->GetLargestPossibleRegion().GetSize()[d] * image->GetSpacing()[d];
    }

    currentTransform->SetCenter(center);


    statismo::ASMPreprocessedImage<MeshType> *pimage = aModel->GetstatismoImplObj()->GetImagePreprocessor()->Preprocess(image);

    // just for testing
    itk::ReducedVarianceModelBuilder<MeshType>::Pointer redModelBuilder = itk::ReducedVarianceModelBuilder<MeshType>::New();

    StatisticalModelType::Pointer reducedMOdel = redModelBuilder->BuildNewModelWithLeadingComponents(aModel->GetStatisticalModel(), 20);
    aModel->SetStatisticalModel(reducedMOdel);
    statismo::VectorType coeffs = statismo::VectorType::Zero(aModel->GetStatisticalModel()->GetNumberOfPrincipalComponents());


    std::vector<PointType> linePoints;
    linePoints = readLandmarksFile<MeshType>(std::string("/tmp/0021lms-line.csv"));

    // Get poitn ids of reference and target points
    std::vector<PointType> refPoints;
    refPoints = readLandmarksFile<MeshType>(std::string("/tmp/fancylms.csv"));

    std::vector<PointType> targetPoints;
    targetPoints = readLandmarksFile<MeshType>(std::string("/tmp/0021lms.csv"));


    typedef itk::PointsLocator< typename RepresenterType::MeshType::PointsContainer > PointsLocatorType;
    PointsLocatorType::Pointer ptLocator = PointsLocatorType::New();
    ptLocator->SetPoints(aModel->GetStatisticalModel()->GetRepresenter()->GetReference()->GetPoints());
    ptLocator->Initialize();

    if (refPoints.size() != targetPoints.size()) {
        std::cout << "need the same number of reference and target points" << std::endl;
        exit(-1);
    }

    FittingStepType::CorrespondencePoints correspondingPoints;
    for (unsigned i = 0; i < refPoints.size(); ++i) {
        long ptId = ptLocator->FindClosestPoint(refPoints[i]);
        correspondingPoints.push_back(std::make_pair(ptId, targetPoints[i]));
    }




    // very ITK unlike, we use a init method instead of setting all fields manually.
    // This avoids 99% of all core dumps :-)
    FittingStepType::Pointer fittingStep = FittingStepType::New();
    fittingStep->init(image, pimage, correspondingPoints, linePoints, aModel, FittingStepType::SamplerPointerType(fitSampler.GetPointer()), mhFitConfig, currentTransform, coeffs);

    std::cout << "Initialization done." << std::endl;


    for (int i =1; i <= 1500; ++i) {

        std::cout << "iteration: " << i << std::endl;

        fittingStep->NextSample();
        FittingResultType::Pointer result = fittingStep->GetOutput();
        if (!result->IsValid()) {
            std::cout << "invalid result, aborting " <<std::endl;
            exit(42);
        }
//        coeffs = fromVnlVector(result->GetCoefficients());
//        std::cout << "coeffs (adj)" << toVnlVector(coeffs) << std::endl;

        std::cout << "rigid params " << result->GetRigidTransformParameters() << std::endl;

        std::cout << "Writing result of iteration " << i << std::endl;
        itk::MeshFileWriter<MeshType>::Pointer writer = itk::MeshFileWriter<MeshType>::New();
        std::stringstream filename;
        filename << "/tmp/itkmesh-" << i << ".vtk";
        writer->SetFileName(filename.str());
        MeshType::Pointer ms = result->GetMesh();
        writer->SetInput(ms);
        writer->Update();

    }




    return 0;
}

