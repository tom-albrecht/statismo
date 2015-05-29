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

#ifndef STATISMO_ASMNORMALDIRECTIONPOINTSAMPLER_H
#define STATISMO_ASMNORMALDIRECTIONPOINTSAMPLER_H

#include <itkObject.h>
#include <itkMacro.h>
#include <itkPointsLocator.h>
#include "itkASMPointSampler.h"
#include "itkTriangleMeshAdapter.h"
#include "itkCovariantVector.h"

namespace itk {
    template <typename TPointSet>
    class ASMNormalDirectionPointSampler: public ASMPointSampler<TPointSet> {

        typedef itk::PointsLocator<typename TPointSet::PointsContainer> PointsLocatorType;
        typedef typename PointsLocatorType::Pointer PointsLocatorPointerType;
        typedef typename PointsLocatorType::PointsContainer PointsContainerType;
        typedef itk::Vector<typename TPointSet::PixelType> VectorType;
        typedef itk::TriangleMeshAdapter<typename TPointSet::PixelType> MeshAdapterType;
        typedef typename MeshAdapterType::Pointer MeshAdapterPointerType;
        typedef typename MeshAdapterType::PointNormalType PointNormalType;
        typedef typename MeshAdapterType::PointNormalsContainerPointer PointNormalsContainerPointer;

    public:
        /* standard ITK typedefs and macros */
        typedef ASMNormalDirectionPointSampler Self;
        typedef ASMPointSampler<TPointSet> Superclass;
        typedef SmartPointer <Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;
        itkNewMacro( Self );
        itkTypeMacro( Self, Object);


        typedef typename statismo::Representer<TPointSet>::PointType PointType;

        unsigned int GetNumberOfPoints() const {
            return m_numberOfPoints;
        }

        void SetNumberOfPoints(unsigned int numberOfPoints) {
            m_numberOfPoints = numberOfPoints;
        }

        float GetPointSpacing() const {
            return m_pointSpacing;
        }

        void SetPointSpacing(float pointSpacing) {
            m_pointSpacing = pointSpacing;
        }

        virtual std::vector<PointType> Sample(const TPointSet* const pointSet, const PointType& samplePoint) const {
            //FIXME: this only works for TriangleMeshes for now.

            //std::cout << "SAMPLING AT " << samplePoint << std::endl;
            std::vector<PointType> samples;
            samples.reserve(m_numberOfPoints);

            unsigned int normalPointId = FindClosestPointId(pointSet, samplePoint);

            PointNormalType normal = GetNormalForPointId(pointSet, normalPointId);
            // Convert to an itk::(ContraVariant)Vector, because no operators are overloaded for adding the
            // covariant vector normal to a point.
            VectorType normalVector(normal.GetDataPointer());

            int startInclusive = -(m_numberOfPoints / 2);
            int endExclusive = (m_numberOfPoints + 1) / 2;

            for(int i = startInclusive; i < endExclusive; ++i) {
                PointType sample = samplePoint + normalVector * i * m_pointSpacing;
                //std::cout << "sample(" << i <<") = " << sample << std::endl;
                samples.push_back(sample);
            }

            return samples;
        }

        const unsigned int FindClosestPointId(const TPointSet* const pointSet, const PointType &targetPoint) const {
            PointsLocatorPointerType locator = PointsLocatorType::New();
            locator->SetPoints(const_cast<PointsContainerType*>(pointSet->GetPoints()));
            locator->Initialize();
            return locator->FindClosestPoint(targetPoint);
        }

        PointNormalType GetNormalForPointId(const TPointSet* const pointSet, const unsigned &targetPointId) const {
            MeshAdapterPointerType adapter = MeshAdapterType::New();
            adapter->SetMesh(const_cast<TPointSet*>(pointSet));
            PointNormalsContainerPointer normals = adapter->GetPointNormals();

            PointNormalType normal = normals->GetElement(targetPointId);
            return normal * (1.0 / normal.GetNorm());
        }

        PointNormalType GetNormalForPoint(const TPointSet* const pointSet, const PointType &targetPoint) const {
            return GetNormalForPointId(pointSet, FindClosestPointId(pointSet, targetPoint));
        }


            private:
        unsigned int m_numberOfPoints;
        float m_pointSpacing;

    };

}
#endif //STATISMO_ASMNORMALDIRECTIONPOINTSAMPLER_H
