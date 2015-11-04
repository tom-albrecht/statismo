/*
 * This file is part of the statismo library.
 *
 * Author: Christoph Langguth (christoph.langguth@unibas.ch)
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
#ifndef MESHREPRESENTER_H_
#define MESHREPRESENTER_H_

#include "Representer.h"

/**
 * \brief Provides the interface between statismo and the dataset type the application uses.
 *
 * A Representer is a type that provides the connection between the statismo library
 * and the application. It distinguishes three different representations of the data, and provides methods for conversion between those representations:
 * - a Dataset, typically as read from a file on the disk
 * - a Sample, which is a geometric (generally a rigid or affine) transform of the dataset
 * - a SampleVector, which is an internal representation (vector) useful from the statistical analysis.
 *
 * In the following the methods and types that have to be implemented to write a new
 * Representer for your application are given.
 *
 * \warning This class is never actually used, but serves only for documentation purposes.
 */
//RB: would it be possible to make all representers inherit from it, so as to strictly enforce the interface?
namespace statismo {

template<class T>
class MeshRepresenter: public Representer<T> {

  public:

    typedef typename RepresenterTraits<T>::RigidTransformType RigidTransformType;
    typedef typename RepresenterTraits<T>::RigidTransformPointerType RigidTransformPointerType;
    typedef typename RepresenterTraits<T>::DatasetPointerType MeshPointerType;

    virtual ~MeshRepresenter() {
    }

    /** Clone the representer */
    virtual MeshRepresenter* Clone() const = 0;

    /** Delete the representer object */
    virtual void Delete() const = 0;

    virtual typename T::PointType TransformPoint(typename T::PointType &point, const RigidTransformPointerType transform, bool inverse = false) const = 0;

    virtual MeshPointerType TransformMesh(MeshPointerType mesh, const RigidTransformPointerType transform) const = 0;

    virtual RigidTransformPointerType ComputeRigidTransformFromLandmarks(const std::vector<typename T::PointType> &fixedLandmarks, const std::vector<typename T::PointType> &movingLandmarks) const = 0;
};

class RigidTransformation {

};

}

#endif /* MESHREPRESENTER_H_ */

