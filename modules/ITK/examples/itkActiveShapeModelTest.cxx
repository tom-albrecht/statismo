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
#include "ASMFeatureExtractor.h"
#include "itkASMNormalDirectionGradientGaussianFeatureExtractor.h"
#include "itkActiveShapeModel.h"
#include "itkASMNormalDirectionPointSampler.h"

typedef itk::Mesh<float, 3> MeshType;
typedef itk::Image<float, 3> ImageType;
typedef itk::ActiveShapeModel<MeshType, ImageType> ActiveShapeModelType;
typedef itk::ImageFileReader<ImageType> ImageReaderType;

typedef itk::StandardMeshRepresenter<float, 3> RepresenterType;
typedef itk::StatisticalModel<MeshType> StatisticalModelType;

int main(int argc, char *argv[]) {
    //statismo::ASMFeatureExtractorFactory<vtkPolyData, vtkStructuredPoints>::RegisterImplementation(vtkASMNormalDirectionFeatureExtractorFactory::GetInstance());
    statismo::ASMFeatureExtractorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMNormalDirectionFeatureExtractorFactory<MeshType, ImageType>::GetInstance());

    //StatisticalModelType::Pointer ssm = StatisticalModelType::New();
    ActiveShapeModelType::Pointer aModel = ActiveShapeModelType::New();
    new ActiveShapeModelType;

    //std::string modelname("/home/langguth/workspaces.exported/stk.idea/bladderdemo/src/main/resources/bladder/asmModels/asmLevel-1.h5");
    std::string modelname("/tmp/asmLevel-1.h5");
    RepresenterType::Pointer representer = RepresenterType::New();
    aModel->Load(representer, modelname.c_str());

    std::cout << "statismo::SSM impl object: " << aModel->GetStatisticalModel()->GetstatismoImplObj() << std::endl;
    std::cout << "statismo::SSM impl object: " << aModel->GetStatisticalModel()->GetstatismoImplObj() << std::endl;
    std::cout << "statismo::SSM impl object: " << aModel->GetStatisticalModel()->GetstatismoImplObj() << std::endl;
    StatisticalModelType::Pointer s = aModel->GetStatisticalModel();

    itk::ASMNormalDirectionPointSampler::Pointer fitSampler = itk::ASMNormalDirectionPointSampler::New();
    fitSampler->SetNumberOfPoints(42);
    fitSampler->SetPointSpacing((float)0.5);

    statismo::ASMFitterConfiguration fitConfig(1,2,3);

//    ImageReaderType::Pointer reader = ImageReaderType::New();
//    reader->SetFileName("/tmp/24.vtk");
//    reader->Update();
//    ImageType::Pointer image = reader->GetOutput();
//    MeshType::Pointer mesh = aModel->GetStatisticalModel()->DrawMean();

    typename itk::ASMFitter::Pointer fitter = itk::ASMFitter::New();
    fitter->SetSourceModel(aModel); // "source"
    fitter->SetTargetImage(image); // "target"
    fitter->SetSamplingStrategy(fitSampler); // "how"
    fitter->SetThresholds(fitConfig); // "success criteria"
    fitter->SetMaxIterations(42); // "how often"

//
//    typename itk::ASMFitter<ASM>::Pointer fitter = itk::ASMFitter<ASM>::New();
//    fitter->Init(fitConfig, aModel, mesh, image, pointSampler);
//    for (int i=0; i < 25; ++i) {
//        std::cout << "iteration " << i << " start" << std::endl;
//        itk::ASMFitterResult::Pointer result = fitter->Fit();
//        if (!result->IsValid()) {
//            std::cout << "invalid result, aborting " <<std::endl;
//            exit(0);
//        }
//        std::cout << "coeffs (adj)" << result->GetCoefficients() << std::endl;
//        fitter = fitter->SetPointset(result->GetMesh());
//    }
    std::cout << std::endl << "#### Program ends here ####" << std::endl << std::endl;
    return 0;
}
//
//
