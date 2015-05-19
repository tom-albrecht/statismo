/*
 * This file is part of the statismo library.
 *
 * Author: Christoph Langguth (christoph.langguth@unibas.ch)
 *
 * Copyright (c) 2011-2015 University of Basel
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

#ifndef STATISMO_ACTIVESHAPEMODEL_H
#define STATISMO_ACTIVESHAPEMODEL_H

#include "MultiVariateNormalDistribution.h"
#include "StatisticalModel.h"
#include "ASMProfile.h"
#include "ASMFeatureExtractor.h"

namespace statismo {

    template <typename TPointSet, typename TImage>
    class ActiveShapeModel {

    typedef Representer<TPointSet> RepresenterType;
    typedef StatisticalModel<TPointSet> StatisticalModelType;
    typedef typename RepresenterType::PointType PointType;
         typedef ASMFeatureExtractorFactory<TPointSet, TImage> FeatureExtractorFactoryType;

        typedef ASMFeatureExtractor<TPointSet, TImage> ASMFeatureExtractorType;

    protected:
        ActiveShapeModel(const StatisticalModelType* statisticalModel, const ASMFeatureExtractorType* fe, std::vector<ASMProfile>& profiles)
        : m_statisticalModel(statisticalModel),
        m_featureExtractor(fe),
        m_profiles(profiles)
        {}

        //virtual StatisticalModelType* LoadStatisticalShapeModel(const RepresenterType* representer, const H5::Group& group) = 0;
    public:

        virtual ~ActiveShapeModel() {
            std::cout << "statismo::ASM destructor" << std:: endl;
        }


        std::vector<ASMProfile> GetProfiles() const {
            return m_profiles;
        }

        const StatisticalModelType* GetStatisticalModel() const {
            return m_statisticalModel;
        }

        const ASMFeatureExtractor<TPointSet, TImage>* GetFeatureExtractor() const {
            return m_featureExtractor;
        }

        statismo::MultiVariateNormalDistribution GetMarginalAtPointId(unsigned pointId) const {
            PointType imean = m_statisticalModel->DrawMeanAtPoint(pointId);
            MatrixType icovariances = m_statisticalModel->GetCovarianceAtPoint(pointId, pointId);

            unsigned int dimensions = m_statisticalModel->GetRepresenter()->GetDimensions();

            statismo::VectorType mean(dimensions);
            statismo::MatrixType covariances(dimensions,dimensions);
            for(int i=0; i< dimensions; ++i) {
                mean[i] = imean[i];
                for(int j=0; j < dimensions; ++j) {
                    covariances(i,j) = icovariances(i,j);
                }
            }

            return statismo::MultiVariateNormalDistribution(mean, covariances);
        }


        static const ActiveShapeModel<TPointSet, TImage>* Load(RepresenterType* representer, const std::string& filename) {
            H5::H5File file;

            try {
                file = H5::H5File(filename, H5F_ACC_RDONLY);
            } catch (H5::Exception &e) {
                std::string msg(std::string("could not open HDF5 file \n") + e.getCDetailMsg());
                throw StatisticalModelException(msg.c_str());
            }

            H5::Group rootGroup = file.openGroup("/");

            StatisticalModelType* statisticalModel = StatisticalModel<TPointSet>::Load(representer, rootGroup);

            H5::Group asmGroup = rootGroup.openGroup("activeShapeModel");
            H5::Group feGroup = asmGroup.openGroup("featureExtractor");
            H5::Group profilesGroup = asmGroup.openGroup("profiles");

            std::vector<int> pointIds;
            statismo::MatrixType means;
            statismo::MatrixType covariances;

            unsigned int numPoints = (unsigned int) statismo::HDF5Utils::readIntAttribute(profilesGroup, "numberOfPoints");
            unsigned int profileLength = (unsigned int) statismo::HDF5Utils::readIntAttribute(profilesGroup, "profileLength");

            statismo::HDF5Utils::readMatrix(profilesGroup, "covariances", covariances);
            statismo::HDF5Utils::readArray(profilesGroup, "pointIds", pointIds);
            statismo::HDF5Utils::readMatrix(profilesGroup, "means", means);

            std::vector<ASMProfile> profiles;
            profiles.reserve(numPoints);

            for (unsigned int i=0; i < numPoints; ++i) {
                unsigned int covOffset = i * profileLength;

                statismo::VectorType mean = means.row(i);

                statismo::MatrixType cov(profileLength, profileLength);
                cov.block(0, 0, profileLength, profileLength) = covariances.block(covOffset, 0, profileLength,
                                                                                  profileLength);


                statismo::MultiVariateNormalDistribution mvd(mean, cov);

                profiles.push_back(statismo::ASMProfile(pointIds[i], mvd));
            }

            std::string feType = HDF5Utils::readStringAttribute(feGroup, "type");
            const FeatureExtractorFactoryType* feFactory = FeatureExtractorFactoryType::GetImplementation(feType);
            if (!feFactory) {
                std::string msg(std::string("No feature extractor implementation found for type: ") + feType);
                throw StatisticalModelException(msg.c_str());
            }
            const ASMFeatureExtractor<TPointSet, TImage>* featureExtractor = feFactory->Instantiate(feGroup);

            ActiveShapeModel* am = new ActiveShapeModel(statisticalModel, featureExtractor, profiles);

            feGroup.close();
            asmGroup.close();
            rootGroup.close();
            file.close();

            return am;
        }

    private:
        std::vector<ASMProfile> m_profiles;
       const StatisticalModelType *m_statisticalModel;
        const ASMFeatureExtractor<TPointSet, TImage>* m_featureExtractor;
    };
};

#endif //STATISMO_ACTIVESHAPEMODEL_H

