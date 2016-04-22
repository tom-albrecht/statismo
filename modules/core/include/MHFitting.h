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


#ifndef STATISMO_MHFITTING_H
#define STATISMO_MHFITTING_H

#include "ActiveShapeModel.h"
#include "ASMFitting.h"
#include "ASMPointSampler.h"
#include <boost/thread.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/thread/future.hpp>

#include <Sampling/DistributionEvaluator.h>
#include <Sampling/MarkovChain.h>
#include <Sampling/Algorithm/ForwardChain.h>
#include <Sampling/Algorithm/Metropolis.h>
#include <Sampling/RandomGenerator.h>
#include <Sampling/ProposalGenerator.h>
#include <Sampling/Proposal/RandomProposal.h>
#include <Sampling/Logger/QuietLogger.h>
#include <Sampling/RandomGenerator.h>
#include <Sampling/Algorithm/MetropolisHastings.h>



namespace statismo {

  template<typename RigidTransformPointerType>
  class MHFittingResult {
  public:
      MHFittingResult() {
      }

      MHFittingResult(VectorType coefficients, RigidTransformPointerType rigidTransform) :
              m_coefficients(coefficients),
              m_rigidTransform(rigidTransform) { }

      VectorType GetCoefficients() {
          return m_coefficients;
      }

      RigidTransformPointerType GetRigidTransform() {
          return m_rigidTransform;
      }

  private:
      VectorType m_coefficients;
      RigidTransformPointerType m_rigidTransform;
  };






  // Proposals
  /*============================================================================*/
  template<class RigidTransformPointerType>
  class GaussianModelUpdate : public sampling::ProposalGenerator<MHFittingResult<RigidTransformPointerType> > {
    private:
      int sigma;

    public:
      GaussianModelUpdate( double stepSize = 0.1 ) : sigma(stepSize) {}

      // ProposalGeneratorInterface interface
    public:
      virtual void generateProposal(MHFittingResult<RigidTransformPointerType>& proposal, const MHFittingResult<RigidTransformPointerType>& currentSample){
        // TODO: Update model parameters with gaussian noise
        proposal = currentSample;
      }
      virtual double transitionProbability(const MHFittingResult<RigidTransformPointerType>& start, const MHFittingResult<RigidTransformPointerType>& end){
        return 0.0;
      }
  };


  // Evaluators
  /*============================================================================*/
  template<class RigidTransformPointerType>
  class PointEvaluator : public sampling::DistributionEvaluator<MHFittingResult<RigidTransformPointerType> > {

      // DistributionEvaluatorInterface interface
    public:
      virtual double evalSample(const MHFittingResult<RigidTransformPointerType>& currentSample){
        return 0.0;
      }
  };


    template<class RigidTransformPointerType>
    class BasicSampling {
      public:
        static sampling::MarkovChain<MHFittingResult<RigidTransformPointerType> >* buildChain(RigidTransformPointerType transform,statismo::VectorType coeffs) {
          // basics
          sampling::RandomGenerator* rGen = new sampling::RandomGenerator(42);
          MHFittingResult<RigidTransformPointerType> init = MHFittingResult<RigidTransformPointerType>(coeffs,transform);

          // Proposals
          GaussianModelUpdate<RigidTransformPointerType>* poseProposalRough = new GaussianModelUpdate<RigidTransformPointerType>(0.1);
          GaussianModelUpdate<RigidTransformPointerType>* poseProposalFine = new GaussianModelUpdate<RigidTransformPointerType>(0.01);
          std::vector< typename sampling::RandomProposal< MHFittingResult<RigidTransformPointerType> >::GeneratorPair> poseProposalsVector(2);
          poseProposalsVector[0] = std::pair<sampling::ProposalGenerator<MHFittingResult<RigidTransformPointerType> >*,double>(poseProposalRough,0.3);
          poseProposalsVector[1] = std::pair<sampling::ProposalGenerator<MHFittingResult<RigidTransformPointerType> >*,double>(poseProposalFine,0.7);
          sampling::RandomProposal<MHFittingResult<RigidTransformPointerType> >* mainProposal = new sampling::RandomProposal<MHFittingResult<RigidTransformPointerType> >(poseProposalsVector,rGen);

          // Evaluators
          PointEvaluator<RigidTransformPointerType>* mainEval = new PointEvaluator<RigidTransformPointerType>();




          sampling::MarkovChain<MHFittingResult<RigidTransformPointerType> >* chain =  new sampling::MetropolisHastings<MHFittingResult<RigidTransformPointerType> >(mainProposal, mainEval, new sampling::QuietLogger<MHFittingResult<RigidTransformPointerType> >(), init, rGen );
          return chain;
        }
    };



    class MHFittingConfiguration {
    private:

    public:
        MHFittingConfiguration(const ASMFittingConfiguration& asmFittingConfiguration) : m_asmFittingConfiguration(asmFittingConfiguration) { }
        const ASMFittingConfiguration& GetAsmFittingconfiguration() const {return m_asmFittingConfiguration; }

    private:
        ASMFittingConfiguration m_asmFittingConfiguration;
    };






    template<typename TPointSet, typename TImage>
    class MHFittingStep {
        typedef ActiveShapeModel<TPointSet, TImage> ActiveShapeModelType;
        typedef ASMPointSampler<TPointSet, TImage> PointSamplerType;
        typedef typename Representer<TPointSet>::PointType PointType;
        typedef StatisticalModel<TPointSet> StatisticalModelType;
        typedef ASMFeatureExtractor<TPointSet, TImage> FeatureExtractorType;
        typedef ASMPreprocessedImage<TPointSet> PreprocessedImageType;
        typedef typename ActiveShapeModelType::RepresenterType::RigidTransformPointerType RigidTransformPointerType;
        typedef  ASMFittingStep<TPointSet, TImage> ASMFittingStepType;

        MHFittingStep(sampling::MarkovChain<MHFittingResult<RigidTransformPointerType> >* chain)
                :
                m_chain(chain) { }

    private:
        class ProfileResult {
        public:
            unsigned int pointId;
            PointType candidatePoint;
            PointType transformedCandidatePoint;
        };


    public:


        static MHFittingStep *Create(sampling::MarkovChain<MHFittingResult<RigidTransformPointerType> >* chain) {
            return new MHFittingStep( chain);
        }

        ~MHFittingStep() {
        }

        MHFittingResult<RigidTransformPointerType> Perform() const {


          //          ASMFittingStepType* asmFittingStep = ASMFittingStepType::Create(m_configuration.GetAsmFittingconfiguration(), m_model, m_sourceCoefficients, m_sourceTransform, m_target, m_sampler);
          //          ASMFittingResult<RigidTransformPointerType> result = asmFittingStep->Perform();


          // runs a markov chain and returns only the accepted proposals
          MHFittingResult<RigidTransformPointerType> params;//(m_sourceTransform,m_sourceCoefficients);
          m_chain->next(params);
          MHFittingResult<RigidTransformPointerType> mhResult(params.GetCoefficients(), params.GetRigidTransform());
          return mhResult;
        }

    private:
        sampling::MarkovChain<MHFittingResult<RigidTransformPointerType> >* m_chain;
    };
}

#endif //STATISMO_ASMFITTING_H

