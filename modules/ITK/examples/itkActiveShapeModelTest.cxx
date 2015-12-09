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


    std::cout << "Initializing..." << std::endl;
    // FIXME: these should go somewhere less "intrusive"
    statismo::ASMFeatureExtractorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMNormalDirectionFeatureExtractorFactory<MeshType, ImageType>::GetInstance());
    statismo::ASMImagePreprocessorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMGaussianGradientImagePreprocessorFactory<MeshType, ImageType>::GetInstance());

    std::string modelname("/export/skulls/data/shapes/submandibular_gland_l/aligned/registered-0248/model-asm/asm-5.h5");
    ActiveShapeModelType::Pointer aModel = ActiveShapeModelType::New();
    RepresenterType::Pointer representer = RepresenterType::New();
    aModel->Load(representer,  modelname.c_str());

    itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::Pointer fitSampler = itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::New();
    fitSampler->SetNumberOfPoints(25);
    fitSampler->SetPointSpacing(1);

    statismo::ASMFittingConfiguration fitConfig(3,5,3);

    ImageReaderType::Pointer reader = ImageReaderType::New();
    //reader->SetFileName("/export/skulls/data/shapes/submandibular_gland_l/aligned/initial/volume-ct/pddca-0522c0002.nii");
    reader->SetFileName("/export/skulls/data/shapes/submandibular_gland_l/raw/normalized/volume-ct/pddca-0522c0002.nii");

    std::vector<PointType> reference;
    std::vector<PointType> target;

    PointType ra,rb,rc,rd,re,rf,ta,tb,tc,td,te,tf;

    ra[0] = 38.3603515625; ra[1] = -61.9775390625; ra[2] = 255.62380981445312; // 38.3603515625,-61.9775390625,255.62380981445312
    rb[0] = 36.21720504760742; rb[1] = -73.27617645263672; rb[2] = 236.9056396484375; // 36.21720504760742,-73.27617645263672,236.9056396484375
    rc[0] = 32.97368621826172; rc[1] = -60.5927734375; rc[2] = 217.5156707763672; // 32.97368621826172,-60.5927734375,217.5156707763672
    rd[0] =  35.740814208984375; rd[1] = -46.073516845703125; rd[2] = 232.21583557128906; // 35.740814208984375,-46.073516845703125,232.21583557128906
    re[0] =  41.47886657714844; re[1] = -60.76667404174805; re[2] = 232.4999237060547; // 41.47886657714844,-60.76667404174805,232.4999237060547
    rf[0] = 29.496631622314453; rf[1] = -61.79425811767578; rf[2] = 235.50245666503906; // 29.496631622314453,-61.79425811767578,235.50245666503906

    reference.push_back(ra);
    reference.push_back(rb);
    reference.push_back(rc);
    reference.push_back(rd);
    reference.push_back(re);
    reference.push_back(rf);

    ta[0] = (float)37.84376907348633; ta[1] = (float)-92.68473052978516; ta[2] = (float)-340.8616638183594; //37.84376907348633,-92.68473052978516,-340.8616638183594
    tb[0] = (float)37.16470718383789; tb[1] = (float)-104.13507843017578; tb[2] = (float)-354.39959716796875; //37.16470718383789,-104.13507843017578,-354.39959716796875
    tc[0] = (float)36.72918701171875; tc[1] = (float)-94.02742004394531; tc[2] = (float)-366.81243896484375; //36.72918701171875,-94.02742004394531,-366.81243896484375
    td[0] = (float)38.515804290771484; td[1] = (float)-78.56915283203125; td[2] = (float)-355.5799255371094; //38.515804290771484,-78.56915283203125,-355.5799255371094
    te[0] = (float)43.56576156616211; te[1] = (float)-93.40294647216797; te[2] = (float)-355.55865478515625; //43.56576156616211,-93.40294647216797,-355.55865478515625
    tf[0] = (float)31.05061912536621; tf[1] = (float)-93.81318664550781; tf[2] = (float)-356.4269714355469; //31.05061912536621,-93.81318664550781,-356.4269714355469

    target.push_back(ta);
    target.push_back(tb);
    target.push_back(tc);
    target.push_back(td);
    target.push_back(te);
    target.push_back(tf);

    RigidTransformType::Pointer currentTransform = 0;
    if (reference.size()) currentTransform = representer->ComputeRigidTransformFromLandmarks(reference, target);

    reader->Update();
    ImageType::Pointer image = reader->GetOutput();
    statismo::ASMPreprocessedImage<MeshType> *pimage = aModel->GetstatismoImplObj()->GetImagePreprocessor()->Preprocess(image);

    // just for testing
    aModel->SetStatisticalModel(aModel->GetStatisticalModel());

    statismo::VectorType coeffs = statismo::VectorType::Zero(aModel->GetStatisticalModel()->GetNumberOfPrincipalComponents());



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

    for (int i =1; i <= 10; ++i) {
        itk::TimeProbe clock;
        clock.Start();
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
        clock.Stop();

        double elapsed = clock.GetMean();
        std::cout << "Elapsed " << elapsed << std::endl;
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

