#include <vtkPolyData.h>
#include <vtkStructuredPoints.h>
#include "ASMFeatureExtractor.h"
#include "ASMPointSampler.h"
#include "ActiveShapeModel.h"
#include "ASMFitting.h"

#include "vtkStandardMeshRepresenter.h"

class vtkASMNormalDirectionFeatureExtractor : public statismo::ASMFeatureExtractor<vtkPolyData, vtkStructuredPoints> {
public:
    const vtkASMNormalDirectionFeatureExtractor* SetImage(const vtkStructuredPoints* image) const {return this;};
    const vtkASMNormalDirectionFeatureExtractor* SetPointset(const vtkPolyData *dataset) const  {return this;};
    statismo::VectorType ExtractFeatures(const statismo::vtkPoint& point) const {return statismo::VectorType::Zero(5);} ;

};

class vtkASMNormalDirectionFeatureExtractorFactory : public statismo::ASMFeatureExtractorFactory<vtkPolyData, vtkStructuredPoints> {
public:

    static const vtkASMNormalDirectionFeatureExtractorFactory* GetInstance() {
        return impl;
    }

    virtual std::string GetDescriptor() const {
        return "vtkASMNormalDirectionFeatureExtractorFactory";
    }

    virtual const statismo::ASMFeatureExtractor<vtkPolyData, vtkStructuredPoints>* Instantiate(
             const H5::Group &h5Group) const {
        return new vtkASMNormalDirectionFeatureExtractor();
    }

private:
    static vtkASMNormalDirectionFeatureExtractorFactory* impl;
};
vtkASMNormalDirectionFeatureExtractorFactory* vtkASMNormalDirectionFeatureExtractorFactory::impl = 0;




class vtkPointSampler : public statismo::ASMPointSampler<vtkPolyData> {
public:
    vtkPointSampler(vtkPolyData* mesh) : ASMPointSampler<vtkPolyData>(mesh) {}

    const vtkPointSampler* SetPointSet(const vtkPolyData* mesh) const {return this;};
    std::vector<statismo::vtkPoint> SampleAtPoint(const statismo::vtkPoint& targetPoint) const {
        std::vector<statismo::vtkPoint> v;
        return v;
    }

    std::vector<statismo::vtkPoint> SampleAtPointId(const unsigned& targetPointId) const {
        std::vector<statismo::vtkPoint> v;
        return v;
    }


private:
    vtkPolyData* m_mesh;
};

int main(int argc, char *argv[]) {

    statismo::ASMFeatureExtractorFactory<vtkPolyData, vtkStructuredPoints>::RegisterImplementation(vtkASMNormalDirectionFeatureExtractorFactory::GetInstance());

    typedef statismo::ActiveShapeModel<vtkPolyData, vtkStructuredPoints> vtkActiveShapeModelType;


    std::string modelname("/home/langguth/workspaces.exported/stk.idea/bladderdemo/src/main/resources/bladder/asmModels/asmLevel-1.h5");
    statismo::vtkStandardMeshRepresenter* representer = statismo::vtkStandardMeshRepresenter::Create();
    vtkStructuredPoints* image = vtkStructuredPoints::New();
    const vtkActiveShapeModelType* model = vtkActiveShapeModelType::Load(representer, modelname);


    vtkPolyData* mesh = model->GetStatisticalModel()->DrawMean();

    vtkPointSampler sampler(mesh);

    statismo::ASMFitterConfiguration fitConfig(5,3, 3);

    statismo::ASMFitter<vtkPolyData, vtkStructuredPoints>* fitter = statismo::ASMFitter<vtkPolyData, vtkStructuredPoints>::Create(fitConfig, model, mesh, image, &sampler);

    for (int i=0; i < 25; ++i) {
        std::cout << "iteration " << i << " start" << std::endl;
        statismo::ASMFitterResult<vtkPolyData, vtkStructuredPoints> result = fitter->Fit();
        if (!result.IsValid()) {
            std::cout << "invalid result, aborting " <<std::endl;
            exit(0);
        }
        std::cout << "coeffs (adj)" << result.GetCoefficients() << std::endl;
        fitter = fitter->SetMesh(result.GetMesh());
    }
    return 0;
}

