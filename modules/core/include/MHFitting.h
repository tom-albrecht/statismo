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
#include <Sampling/Evaluator/ProductEvaluator.h>
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
  using namespace sampling;
  using namespace std;

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


  template< class TPointSet, class TImage>
  class mcmc {


      typedef ActiveShapeModel<TPointSet, TImage> ActiveShapeModelType;
      typedef ASMPointSampler<TPointSet, TImage> PointSamplerType;
      typedef typename Representer<TPointSet>::PointType PointType;
      typedef StatisticalModel<TPointSet> StatisticalModelType;
      typedef typename StatisticalModelType::DatasetPointerType SampleType;
      typedef typename StatisticalModelType::ValueType ValueType;
      typedef ASMFeatureExtractor<TPointSet, TImage> FeatureExtractorType;
      typedef ASMPreprocessedImage<TPointSet> PreprocessedImageType;
      typedef typename ActiveShapeModelType::RepresenterType::RigidTransformPointerType RigidTransformPointerType;
      typedef  ASMFittingStep<TPointSet, TImage> ASMFittingStepType;
      typedef MHFittingResult<RigidTransformPointerType> ChainSampleType;


      // Proposals
      /*============================================================================*/
      class GaussianModelUpdate : public ProposalGenerator<ChainSampleType > {
        private:
          int sigma;

        public:
          GaussianModelUpdate( double stepSize = 0.1 ) : sigma(stepSize) {}

          // ProposalGeneratorInterface interface
        public:
          virtual void generateProposal(ChainSampleType& proposal, const ChainSampleType& currentSample){
            VectorType shapeParams = currentSample.GetCoefficients();

            // TODO: Update model parameters with gaussian noise
            proposal = ChainSampleType(shapeParams,currentSample.GetRigidTransform());
          }
          virtual double transitionProbability(const ChainSampleType& start, const ChainSampleType& end){
            return 0.0;
          }
      };


      class ASMModelUpdate : public ProposalGenerator<ChainSampleType > {
        private:
          int N;
          ActiveShapeModelType* asmodel;
          ASMFittingConfiguration* fittingConfiguration;
          TImage* target;
          PointSamplerType* sampler;

        public:
          ASMModelUpdate( ActiveShapeModelType* asmodel, int nSteps = 1 ) : N(nSteps), asmodel(asmodel) {}

          // ProposalGeneratorInterface interface
        public:
          virtual void generateProposal(ChainSampleType& proposal, const ChainSampleType& currentSample){

            ASMFittingStepType* asmFittingStep = ASMFittingStepType::Create(fittingConfiguration, asmodel, currentSample.GetCoefficients(), currentSample.GetRigidTransform(), target, sampler);
            ASMFittingResult<RigidTransformPointerType> result = asmFittingStep->Perform();
            proposal = ChainSampleType(result.GetCoefficients(),result.GetRigidTransform());
          }
          virtual double transitionProbability(const ChainSampleType& start, const ChainSampleType& end){
            return 0.0;
          }
      };


      // Evaluators
      /*============================================================================*/
      class PositionEvaluator : public DistributionEvaluator<PointType>
      {
        public:
          /// dummy gradient, zero, do nothing
          virtual void evalGradient( PointType& gradient, PointType const& sample )
          {
            throw logic_error("no gradient implemented in PositionEvaluator");
          }

          /// Returns whether this evaluator is separable into x and y coordinates
          virtual bool isSeparable() const
          {
            return false;
          }
      };

      class Gaussian3DPositionDifferenceEvaluator : public DistributionEvaluator<PointType>
      {
        public:
          Gaussian3DPositionDifferenceEvaluator( double sigma = 1.0) {
            sigma = sigma;
            dNormalizer = -1.5*log( 2*M_PI ) - 3.0 * log( sigma );
          }


          double evalSample( PointType const& diff )
          {
            double dVal = diff.norm2; // TODO check if PointType has this function norm2
            return -dVal/(2.0*sigma*sigma) + dNormalizer;
          }

          void evalGradient( PointType& grad, PointType const& diff )
          {
            grad = -diff / (sigma * sigma);
          }

          virtual bool isSeparable() const
          {
            return true;
          }

        private:
          double sigma;
          double dNormalizer;
      };


      class PointEvaluator : public DistributionEvaluator< ChainSampleType > {
        public:
          PointEvaluator( vector< pair<PointType, int> > targetPointsWithModelIndices, ActiveShapeModelType* asmodel, PositionEvaluator* evaluator) :
            pointsWithIndex(targetPointsWithModelIndices),
            smodel(asmodel->GetStatisticalModel()),
            eval(evaluator)
          {}

          // DistributionEvaluatorInterface interface
        public:
          virtual double evalSample(const ChainSampleType& currentSample){
            double val = 0.0;

            // TODO to draw full sample only for landmarks is inefficient
            SampleType sample = smodel->DrawSample(currentSample.GetCoefficients());
            for( int i = 0; i < pointsWithIndex.size(); ++i) {
              PointType genPoint = smodel->GetRepresenter()->PointSampleFromSample(sample,pointsWithIndex(i).second);
              val += eval->evalSample(genPoint-pointsWithIndex[i].first); // TODO: check weather we can subtract two points from each others and get a point back
            }
            return 0.0;
          }
        private:
          vector< pair<PointType, int> > pointsWithIndex;
          StatisticalModelType* smodel;
          PositionEvaluator* eval;

      };


      // TODO: Full image evaluator, i.e. ASM Evaluator is missing


    public:
      // "Script"
      /*============================================================================*/
      class BasicSampling {
        public:
          static MarkovChain<ChainSampleType >* buildChain(
              vector<PointType,int> targetPointsWithIndex,
              ActiveShapeModelType* asmodel,
              RigidTransformPointerType transform,
              statismo::VectorType coeffs) {

            // basics
            RandomGenerator* rGen = new RandomGenerator(42);
            ChainSampleType init = ChainSampleType(coeffs,transform);

            // Proposals // TODO: add ASMProposal
            GaussianModelUpdate* poseProposalRough = new GaussianModelUpdate(0.1);
            GaussianModelUpdate* poseProposalFine = new GaussianModelUpdate(0.01);
            vector< typename RandomProposal< ChainSampleType >::GeneratorPair> poseProposalsVector(2);
            poseProposalsVector[0] = pair<ProposalGenerator<ChainSampleType >*,double>(poseProposalRough,0.3);
            poseProposalsVector[1] = pair<ProposalGenerator<ChainSampleType >*,double>(poseProposalFine,0.7);
            RandomProposal<ChainSampleType >* mainProposal = new RandomProposal<ChainSampleType >(poseProposalsVector,rGen);

            // Evaluators // TODO: integrate ASM evaluator
            Gaussian3DPositionDifferenceEvaluator* diffEval = new Gaussian3DPositionDifferenceEvaluator(0.1);
            PointEvaluator* pointEval = new PointEvaluator(targetPointsWithIndex,asmodel,diffEval);

            Gaussian3DPositionDifferenceEvaluator* wideDiffEval = new Gaussian3DPositionDifferenceEvaluator(0.5);
            PointEvaluator* widePointEval = new PointEvaluator(targetPointsWithIndex,asmodel,wideDiffEval);

            std::vector<DistributionEvaluator<ChainSampleType >*> evaluatorList(2);
            evaluatorList[0] = pointEval;
            evaluatorList[1] = widePointEval;
            ProductEvaluator<ChainSampleType >* productEvaluator = new ProductEvaluator<ChainSampleType >(evaluatorList);

            // markov chain
            MarkovChain<ChainSampleType >* chain =  new MetropolisHastings<ChainSampleType >(mainProposal, productEvaluator, new QuietLogger<ChainSampleType >(), init, rGen );
            return chain;
          }
      };

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
      typedef MHFittingResult<RigidTransformPointerType> ChainSampleType;

      MHFittingStep(MarkovChain<ChainSampleType >* chain)
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


      static MHFittingStep *Create(MarkovChain<ChainSampleType >* chain) {
        return new MHFittingStep( chain);
      }

      ~MHFittingStep() {
      }

      ChainSampleType Perform() const {


        //          ASMFittingStepType* asmFittingStep = ASMFittingStepType::Create(m_configuration.GetAsmFittingconfiguration(), m_model, m_sourceCoefficients, m_sourceTransform, m_target, m_sampler);
        //          ASMFittingResult<RigidTransformPointerType> result = asmFittingStep->Perform();


        // runs a markov chain and returns only the accepted proposals
        ChainSampleType params;//(m_sourceTransform,m_sourceCoefficients);
        m_chain->next(params);
        ChainSampleType mhResult(params.GetCoefficients(), params.GetRigidTransform());
        return mhResult;
      }

    private:
      MarkovChain<ChainSampleType >* m_chain;
  };
}

#endif //STATISMO_ASMFITTING_H

