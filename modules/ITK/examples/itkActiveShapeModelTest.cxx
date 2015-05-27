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
#include "ASMFeatureExtractor.h"
#include "itkASMNormalDirectionGradientGaussianFeatureExtractor.h"
#include "itkActiveShapeModel.h"

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

    std::string modelname("/home/langguth/workspaces.exported/stk.idea/bladderdemo/src/main/resources/bladder/asmModels/asmLevel-1.h5");
    RepresenterType::Pointer representer = RepresenterType::New();
    aModel->Load(representer, modelname.c_str());

    std::cout << "statismo::SSM impl object: " << aModel->GetStatisticalModel()->GetstatismoImplObj() << std::endl;
    std::cout << "statismo::SSM impl object: " << aModel->GetStatisticalModel()->GetstatismoImplObj() << std::endl;
    std::cout << "statismo::SSM impl object: " << aModel->GetStatisticalModel()->GetstatismoImplObj() << std::endl;
    StatisticalModelType::Pointer s = aModel->GetStatisticalModel();

    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName("/home/langguth/workspaces.exported/stk.idea/bladderdemo/src/main/resources/bladder/volumes/24.vtk");
    reader->Update();
    ImageType::Pointer image = reader->GetOutput();

//
//
    MeshType::Pointer mesh = aModel->GetStatisticalModel()->DrawMean();
//
//    itk::ASMNormalDirectionPointSampler<ASM>::Pointer sampler = itk::ASMNormalDirectionPointSampler<ASM>::New();
//    //sampler->Init(90, (float)(2.0/3.0));
//    sampler->Init(25, (float)(0.5));
//    //sampler->Init(125, (float)(0.5));
//
//    itk::ASMFitterConfiguration::Pointer fitConfig = itk::ASMFitterConfiguration::New();
//    fitConfig->Set1()
//    fitConfig->Set1()
//    fitConfig->Set1()

//    itk::Point
//
//    typename itk::ASMFitter<ASM>::Pointer fitter = itk::ASMFitter<ASM>::New();
//    typename itk::ASMNormalDirectionPointSampler<ASM>::Pointer pointSampler = itk::ASMNormalDirectionPointSampler<ASM>::New();
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
    std::cout << "Program ends" << std::endl;
    return 0;
}
//
//
