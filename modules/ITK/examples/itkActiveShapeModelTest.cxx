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
#include "itkASMNormalDirectionFeatureExtractor.h"
#include "itkActiveShapeModel.h"
#include "itkASMNormalDirectionPointSampler.h"
#include "itkASMFitter.h"

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
    std::string modelname("/tmp/newasm/asmLevel-1.h5");
    RepresenterType::Pointer representer = RepresenterType::New();
    aModel->Load(representer, modelname.c_str());

    StatisticalModelType::Pointer ssm = aModel->GetStatisticalModel();

    itk::ASMNormalDirectionPointSampler<MeshType>::Pointer fitSampler = itk::ASMNormalDirectionPointSampler<MeshType>::New();
    fitSampler->SetNumberOfPoints(3);
    fitSampler->SetPointSpacing(1);

    statismo::ASMFitterConfiguration fitConfig(30,30,40);

    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName("/tmp/24.vtk");
    reader->Update();
    ImageType::Pointer image = reader->GetOutput();

    //FIXME: yuck
    MeshType::Pointer mean = ssm->DrawMean();
    StatisticalModelType::VectorType mmean = ssm->ComputeCoefficientsForDataset(mean);
    statismo::VectorType smean = Eigen::Map<const statismo::VectorType>(mmean.data_block(), mmean.size());

    typedef itk::ASMFitterStep<MeshType, ImageType> FitterStepType;
    typedef typename FitterStepType::Pointer FitterStepPointerType;
    FitterStepPointerType fitter = FitterStepType::New();
    fitter->SetModel(aModel); // "source"
    fitter->SetSource(smean); // "source"
    fitter->SetTarget(image); // "target"
    fitter->SetSampler(FitterStepType::SamplerPointerType(fitSampler.GetPointer())); // "how"
    fitter->SetConfiguration(fitConfig); // "success criteria"
    fitter->Update();
    statismo::ASMFitterResult result = fitter->GetOutput();

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
