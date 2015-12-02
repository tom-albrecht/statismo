/*
 * This file is part of the statismo library.
 *
 * Author: Marcel Luethi (marcel.luethi@unibas.ch)
 *
 * Copyright (c) 2011 University of Basel
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the project's author nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef ITKSTATISTICALMODEL_H_
#define ITKSTATISTICALMODEL_H_

#include <boost/bind.hpp>
#include <boost/utility/result_of.hpp>

#include <itkObject.h>
#include <itkObjectFactory.h>
#include <itkVersorRigid3DTransform.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "ModelInfo.h"
#include "Representer.h"
#include "StatisticalModel.h"
#include "statismoITKConfig.h"

namespace itk {

    /**
 * \brief ITK Wrapper for the statismo::StatisticalModel class.
 * \see statismo::StatisticalModel for detailed documentation.
 */
template <class T>
class StatisticalModel : public Object {
  public:

    typedef StatisticalModel            Self;
    typedef Object	Superclass;
    typedef SmartPointer<Self>                Pointer;
    typedef SmartPointer<const Self>          ConstPointer;

    itkNewMacro( Self );
    itkTypeMacro( StatisticalModel, Object );

    typedef statismo::Representer<T> RepresenterType;

    // statismo stuff
    typedef statismo::StatisticalModel<T> ImplType;

    typedef typename statismo::DataManager<T>::DataItemType     DataItemType;

    typedef vnl_matrix<statismo::ScalarType> MatrixType;
    typedef vnl_vector<statismo::ScalarType> VectorType;

    template <class F>
    typename boost::result_of<F()>::type callstatismoImpl(F f) const {
        try {
            return f();
        } catch (statismo::StatisticalModelException& s) {
            itkExceptionMacro(<< s.what());
        }
    }

    virtual void SetstatismoImplObj(ImplType* impl) {
        if (m_impl) {
            delete m_impl;
        }
        m_impl = impl;
    }

    virtual ImplType* GetstatismoImplObj() const {
        return m_impl;
    }

    StatisticalModel() : m_impl(0) {}

    virtual ~StatisticalModel() {
        if (m_impl) {
            delete m_impl;
        }
    }


    typedef typename RepresenterType::DatasetPointerType DatasetPointerType;
    typedef typename RepresenterType::DatasetConstPointerType DatasetConstPointerType;

    typedef typename RepresenterType::ValueType ValueType;
    typedef typename RepresenterType::PointType PointType;

    typedef typename statismo::StatisticalModel<T>::PointValuePairType PointValuePairType;
    typedef typename statismo::StatisticalModel<T>::PointValueListType PointValueListType;

    typedef typename statismo::StatisticalModel<T>::PointCovarianceMatrixType PointCovarianceMatrixType;
    typedef typename statismo::StatisticalModel<T>::PointValueWithCovariancePairType PointValueWithCovariancePairType;
    typedef typename statismo::StatisticalModel<T>::PointValueWithCovarianceListType PointValueWithCovarianceListType;

    typedef typename statismo::StatisticalModel<T>::DomainType DomainType;

    void Load(RepresenterType* representer, const char* filename) {
        try {
            SetstatismoImplObj(ImplType::Load(representer, filename));
        } catch (statismo::StatisticalModelException& s) {
            itkExceptionMacro(<< s.what());
        }
    }


    void Load(RepresenterType* representer, const H5::Group& modelRoot) {
        try {
            SetstatismoImplObj(ImplType::Load(representer, modelRoot));
        } catch (statismo::StatisticalModelException& s) {
            itkExceptionMacro(<< s.what());
        }
    }

//    Pointer Transform(RigidTransformPointerType rigidTransform) const {
//        Pointer p = Self::New();
//        BasicRigidTransformationType* basicTransform = new BasicRigidTransformationType(rigidTransform);
//        p->SetstatismoImplObj(m_impl->Transform(basicTransform));
//        delete basicTransform;
//        return p;
//    }


        //TODO: wrap StatisticalModel* BuildReducedVarianceModel( double pcvar );

    const RepresenterType* GetRepresenter() const {
        return callstatismoImpl(boost::bind(&ImplType::GetRepresenter, this->GetstatismoImplObj()));
    }

    const DomainType& GetDomain() const {
        return callstatismoImpl(boost::bind(&ImplType::GetDomain, this->GetstatismoImplObj()));
    }

    DatasetPointerType DrawMean() const {
        return callstatismoImpl(boost::bind(&ImplType::DrawMean, this->GetstatismoImplObj()));
    }

    ValueType DrawMeanAtPoint(const PointType& pt) const {
        typedef ValueType (ImplType::*functype)(const PointType&) const;
        return callstatismoImpl(boost::bind(static_cast<functype>(&ImplType::DrawMeanAtPoint), this->GetstatismoImplObj(), pt));
    }

    ValueType DrawMeanAtPoint(unsigned ptid) const {
        typedef ValueType (ImplType::*functype)(unsigned) const;
        return callstatismoImpl(boost::bind(static_cast<functype>(&ImplType::DrawMeanAtPoint), this->GetstatismoImplObj(), ptid));
    }

    DatasetPointerType DrawSample(const VectorType& coeffs, bool addNoise = false) const {
        typedef DatasetPointerType (ImplType::*functype)(const statismo::VectorType&, bool) const;
        return callstatismoImpl(boost::bind(static_cast<functype>(&ImplType::DrawSample), this->GetstatismoImplObj(), fromVnlVector(coeffs), addNoise));
    }

    DatasetPointerType DrawSample(bool addNoise = false) const {
        typedef DatasetPointerType (ImplType::*functype)(bool) const;
        return callstatismoImpl(boost::bind(static_cast<functype>(&ImplType::DrawSample), this->GetstatismoImplObj(), addNoise));
    }

    DatasetPointerType DrawPCABasisSample(unsigned componentNumber) const {
        typedef DatasetPointerType (ImplType::*functype)(unsigned) const;
        return callstatismoImpl(boost::bind(static_cast<functype>(&ImplType::DrawPCABasisSample), this->GetstatismoImplObj(), componentNumber));
    }

    ValueType DrawSampleAtPoint(const VectorType& coeffs, const PointType& pt, bool addNoise = false) const {
        typedef ValueType (ImplType::*functype)(const statismo::VectorType&, const PointType&, bool) const;
        return callstatismoImpl(boost::bind(static_cast<functype>(&ImplType::DrawSampleAtPoint), this->GetstatismoImplObj(), fromVnlVector(coeffs), pt, addNoise));
    }

    ValueType DrawSampleAtPoint(const VectorType& coeffs, unsigned ptid, bool addNoise  = false) const  {
        typedef ValueType (ImplType::*functype)(const statismo::VectorType&, unsigned, bool) const;
        return callstatismoImpl(boost::bind(static_cast<functype>(&ImplType::DrawSampleAtPoint), this->GetstatismoImplObj(), fromVnlVector(coeffs), ptid, addNoise));
    }


    VectorType ComputeCoefficientsForDataset(DatasetConstPointerType ds) const {
        return toVnlVector(callstatismoImpl(boost::bind(&ImplType::ComputeCoefficientsForDataset, this->GetstatismoImplObj(), ds)));
    }

    VectorType ComputeCoefficientsForSample(DatasetConstPointerType ds) const {
        return toVnlVector(callstatismoImpl(boost::bind(&ImplType::ComputeCoefficientsForSample, this->GetstatismoImplObj(), ds)));
    }

    VectorType ComputeCoefficientsForDataSample(const DataItemType* sample) const {
        return toVnlVector(callstatismoImpl(boost::bind(&ImplType::ComputeCoefficientsForDataSample, this->GetstatismoImplObj(), sample)));
    }

    double ComputeLogProbabilityOfDataset(DatasetConstPointerType ds) const {
        return callstatismoImpl(boost::bind(&ImplType::ComputeLogProbabilityOfDataset, this->GetstatismoImplObj(), ds));
    }

    double ComputeProbabilityOfDataset(DatasetConstPointerType ds) const {
        return callstatismoImpl(boost::bind(&ImplType::ComputeProbabilityOfDataset, this->GetstatismoImplObj(), ds));
    }

    double ComputeLogProbabilityOfCoefficients(const VectorType& coeffs) const {
        return callstatismoImpl(boost::bind(&ImplType::ComputeLogProbabilityOfCoefficients, this->GetstatismoImplObj(), fromVnlVector(coeffs)));
    }

    double ComputeProbabilityOfCoefficients(const VectorType& coeffs) const {
        return callstatismoImpl(boost::bind(&ImplType::ComputeProbabilityOfCoefficients, this->GetstatismoImplObj(), fromVnlVector(coeffs)));
    }

    double ComputeMahalanobisDistanceForDataset(DatasetConstPointerType ds) const {
        return callstatismoImpl(boost::bind(&ImplType::ComputeMahalanobisDistanceForDataset, this->GetstatismoImplObj(), ds));
    }

    VectorType ComputeCoefficientsForPointValues(const PointValueListType& pvlist, double variance) const {
        typedef statismo::VectorType (ImplType::*functype)(const PointValueListType&, double) const;
        return toVnlVector(callstatismoImpl(boost::bind(static_cast<functype>(&ImplType::ComputeCoefficientsForPointValues), this->GetstatismoImplObj(), pvlist, variance)));
    }

    VectorType ComputeCoefficientsForPointValuesWithCovariance(const PointValueWithCovarianceListType& pvclist) const {
      typedef statismo::VectorType(ImplType::*functype)(const PointValueWithCovarianceListType&) const;
      return toVnlVector(callstatismoImpl(boost::bind(static_cast<functype>(&ImplType::ComputeCoefficientsForPointValuesWithCovariance), this->GetstatismoImplObj(), pvclist)));
    }

    DatasetPointerType DatasetToSample(DatasetConstPointerType ds) const {
        return callstatismoImpl(boost::bind(&ImplType::DatasetToSample, this->GetstatismoImplObj(), ds));
    }

    unsigned GetNumberOfPrincipalComponents() const {
        return callstatismoImpl(boost::bind(&ImplType::GetNumberOfPrincipalComponents, this->GetstatismoImplObj()));
    }

    void Save(const char* modelname) {
        typedef void (ImplType::*functype)(const std::string&) const;
        callstatismoImpl(boost::bind(static_cast<functype>(&ImplType::Save), this->GetstatismoImplObj(), modelname));
    }

    void Save(const H5::Group& modelRoot) {
        typedef void (ImplType::*functype)(const H5::Group&) const;
        callstatismoImpl(boost::bind(static_cast<functype>(&ImplType::Save), this->GetstatismoImplObj(), modelRoot));
    }

    float GetNoiseVariance() const {
        return callstatismoImpl(boost::bind(&ImplType::GetNoiseVariance, this->GetstatismoImplObj()));
    }

    MatrixType GetCovarianceAtPoint(const PointType& pt1, const PointType& pt2) const {
        typedef statismo::MatrixType (ImplType::*functype)(const PointType&, const PointType&) const;
        return  toVnlMatrix(callstatismoImpl(boost::bind(static_cast<functype>(&ImplType::GetCovarianceAtPoint), this->GetstatismoImplObj(), pt1, pt2)));
    }

    MatrixType GetCovarianceAtPoint(unsigned ptid1, unsigned  ptid2) const {
        typedef statismo::MatrixType (ImplType::*functype)(unsigned, unsigned ) const;
        return toVnlMatrix(callstatismoImpl(boost::bind(static_cast<functype>(&ImplType::GetCovarianceAtPoint),this->GetstatismoImplObj(), ptid1, ptid2)));
    }

    MatrixType GetJacobian(const PointType& pt) const {
        typedef statismo::MatrixType (ImplType::*functype)(const PointType&) const;
        return toVnlMatrix(callstatismoImpl(boost::bind(static_cast<functype>(&ImplType::GetJacobian), this->GetstatismoImplObj(), pt)));
    }

    MatrixType GetJacobian(unsigned ptId) const {
        typedef statismo::MatrixType (ImplType::*functype)(unsigned) const;
        return toVnlMatrix(callstatismoImpl(boost::bind(static_cast<functype>(&ImplType::GetJacobian), this->GetstatismoImplObj(), ptId)));
    }

    MatrixType GetPCABasisMatrix() const {
        return toVnlMatrix(callstatismoImpl(boost::bind(&ImplType::GetPCABasisMatrix, this->GetstatismoImplObj())));
    }

    MatrixType GetOrthonormalPCABasisMatrix() const {
        return toVnlMatrix(callstatismoImpl(boost::bind(&ImplType::GetOrthonormalPCABasisMatrix, this->GetstatismoImplObj())));
    }

    VectorType GetPCAVarianceVector() const {
        return toVnlVector(callstatismoImpl(boost::bind(&ImplType::GetPCAVarianceVector, this->GetstatismoImplObj())));
    }

    VectorType GetMeanVector() const {
        return toVnlVector(callstatismoImpl(boost::bind(&ImplType::GetMeanVector, this->GetstatismoImplObj())));
    }

    const statismo::ModelInfo& GetModelInfo() const {
        return callstatismoImpl(boost::bind(&ImplType::GetModelInfo, this->GetstatismoImplObj()));
    }

  private:

    static MatrixType toVnlMatrix(const statismo::MatrixType& M) {
        return MatrixType(M.data(), M.rows(), M.cols());

    }

    static VectorType toVnlVector(const statismo::VectorType& v) {
        return VectorType(v.data(), v.rows());

    }

    static statismo::VectorType fromVnlVector(const VectorType& v) {
        return Eigen::Map<const statismo::VectorType>(v.data_block(), v.size());

    }

    StatisticalModel(const StatisticalModel& orig);
    StatisticalModel& operator=(const StatisticalModel& rhs);

    ImplType* m_impl;
};

}

#endif /* ITKSTATISTICALMODEL_H_ */
