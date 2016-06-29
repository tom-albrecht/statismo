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
#include <Sampling/Logger/BestMatchLogger.h>
#include <Sampling/Proposal/MarkovChainProposal.h>


namespace statismo {
  using namespace sampling;
  using namespace std;





    class MHFittingConfiguration {
    private:

    public:
        MHFittingConfiguration(const ASMFittingConfiguration& asmFittingConfiguration, unsigned numberOfProfilePoints = 100 ) : m_asmFittingConfiguration(asmFittingConfiguration) { }
        const ASMFittingConfiguration& GetAsmFittingconfiguration() const {return m_asmFittingConfiguration; }

    private:
        ASMFittingConfiguration m_asmFittingConfiguration;
    };







  class MHFittingParameters {
    public:
      MHFittingParameters() {
      }

      MHFittingParameters(VectorType coefficients, VectorType rigidParameters) :
        m_coefficients(coefficients),
        m_rigidParameters(rigidParameters) { }

      VectorType GetCoefficients() const {
        return m_coefficients;
      }

      VectorType GetRigidTransformParameters() const {
        return m_rigidParameters;
      }


      int size() {
        // TODO: if this is really used we should change to a different vectorized representation of all parameters
        return m_coefficients.size() + m_rigidParameters.size();
      }

      float& operator[]( int i ) {
        // TODO: if this is really used we should change to a different vectorized representation of all parameters
        throw new std::runtime_error("operator[] not implemented for MHFittingResult");
        return m_coefficients(0,0);
      }

    private:
      VectorType m_coefficients;
        VectorType m_rigidParameters;
  };



    template <class MeshType, class PointType>
    class MeshOperations {
    public:
        virtual std::pair<PointType, long> findClosestPoint(MeshType mesh, PointType pt) const = 0;
        virtual VectorType normalAtPoint(MeshType mesh, PointType pt) const = 0;
        virtual short huAtPoint(MeshType mesh, long ptId) const = 0;
        virtual MeshType transformMesh(const MHFittingParameters& fittingParameters) const = 0;
    };



  template< class TPointSet, class TImage>
  class mcmc {


      typedef ActiveShapeModel<TPointSet, TImage> ActiveShapeModelType;
      typedef ASMPointSampler<TPointSet, TImage> PointSamplerType;
      typedef typename Representer<TPointSet>::PointType PointType;
      typedef StatisticalModel<TPointSet> StatisticalModelType;
      typedef typename StatisticalModelType::RepresenterType RepresenterType;
      typedef typename StatisticalModelType::DatasetPointerType SampleType;
      typedef typename StatisticalModelType::ValueType ValueType;
      typedef ASMFeatureExtractor<TPointSet, TImage> FeatureExtractorType;
      typedef ASMPreprocessedImage<TPointSet> PreprocessedImageType;
      typedef typename ActiveShapeModelType::RepresenterType::RigidTransformPointerType RigidTransformPointerType;
      typedef  ASMFittingStep<TPointSet, TImage> ASMFittingStepType;


      // Proposals
      /*============================================================================*/
      class GaussianModelUpdate : public ProposalGenerator<MHFittingParameters > {
        private:
          double sigmaShape;
          double sigmaPose;
          RandomGenerator* rgen;

        public:
          GaussianModelUpdate( double stepSizeShape, double stepSizePose, RandomGenerator* rgen ) : sigmaShape(stepSizeShape), sigmaPose(stepSizePose), rgen(rgen) {}

          // ProposalGeneratorInterface interface
        public:
          virtual void generateProposal(MHFittingParameters& proposal, const MHFittingParameters& currentSample){
            VectorType shapeParams = currentSample.GetCoefficients();

              VectorType newShapeParams(shapeParams.size());
            for (unsigned i = 0; i < shapeParams.size(); ++i) {
                newShapeParams[i] = shapeParams[i] + rgen->normalDbl() * sigmaShape;
            }


              VectorType newRigidParams = currentSample.GetRigidTransformParameters();
              for (unsigned i = 0; i < currentSample.GetRigidTransformParameters().size(); ++i) {
                  if (i < 3) {// rotation parameters
                    newRigidParams[i] += rgen->normalDbl() * 0.01 * sigmaPose;
                  }
                  else { // tranlation scale
                      newRigidParams[i] += rgen->normalDbl() * sigmaPose;
                  }
              }
            proposal = MHFittingParameters(newShapeParams, newRigidParams);
          }

          virtual double transitionProbability(const MHFittingParameters& start, const MHFittingParameters& end){
            return 0.0;
          }
      };


      /*
      class ASMModelUpdate : public ProposalGenerator<MHFittingParameters > {
        private:
          int N;
          ActiveShapeModelType* m_asmodel;
          ASMFittingConfiguration m_fittingConfiguration;
          PreprocessedImageType* m_target;
          PointSamplerType* m_sampler;

        public:
          ASMModelUpdate(ASMFittingConfiguration config, ActiveShapeModelType* asmodel, PreprocessedImageType* target,  PointSamplerType* sampler, int nSteps = 1 )
                  : m_fittingConfiguration(config),
                  m_asmodel(asmodel) ,
                    m_target(target),
                    m_sampler(sampler),
                  N(nSteps) {}

          // ProposalGeneratorInterface interface
        public:
          virtual void generateProposal(MHFittingParameters& proposal, const MHFittingParameters& currentSample){

            ASMFittingStepType* asmFittingStep = ASMFittingStepType::Create(m_fittingConfiguration, m_asmodel, currentSample.GetCoefficients(), currentSample.GetRigidTransform(), m_target, m_sampler);
            ASMFittingResult<RigidTransformPointerType> result = asmFittingStep->Perform();
              statismo::VectorType newCoeffs = result.GetCoefficients();
              statismo::VectorType currCoeffs = currentSample.GetCoefficients();
              statismo::VectorType newProposal = currCoeffs + (newCoeffs - currCoeffs) * 0.1;
            proposal = MHFittingParameters(newProposal,result.GetRigidTransform());
          }
          virtual double transitionProbability(const MHFittingParameters& start, const MHFittingParameters& end){
            return 0.0;
          }
      };
       */

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


      class ModelPriorEvaluator : public DistributionEvaluator<MHFittingParameters>
      {
      public:
          virtual double evalSample(const MHFittingParameters& currentSample) {

                VectorType coefficients = currentSample.GetCoefficients();
              return -0.5 * coefficients.size() *log( 2*M_PI ) - 0.5 * coefficients.squaredNorm();

          }

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


      template <class T>
      class InLungOrBoneEvaluator : public DistributionEvaluator< MHFittingParameters > {
          typedef MeshOperations<typename Representer<T>::DatasetPointerType, typename Representer<T>::PointType> ClosestPointType;

      public:
          InLungOrBoneEvaluator( const Representer<T>* representer, const ClosestPointType* closestPoint, ActiveShapeModelType* asmodel) :
                  m_asmodel(asmodel),
                  m_closestPoint(closestPoint)

          {}

          // DistributionEvaluatorInterface interface
      public:
          virtual double evalSample(const MHFittingParameters& currentSample) {
              double distance = 0.0;

              // TODO to draw full sample only for landmarks is inefficient
//              typename RepresenterType::DatasetPointerType  sampleShape = m_asmodel->GetStatisticalModel()->DrawSample(currentSample.GetCoefficients());
//              typename RepresenterType::DatasetPointerType sample = m_asmodel->GetRepresenter()->TransformMesh(sampleShape, currentSample.GetRigidTransform());


              //std::cout << currentSample.GetRigidTransform()->GetParameters() << std::endl;

              typename RepresenterType::DatasetPointerType sample = m_closestPoint->TransformMesh(m_asmodel->GetStatisticalModel(), currentSample);
              typename RepresenterType::DomainType::DomainPointsListType  points = m_asmodel->GetRepresenter()->GetDomain().GetDomainPoints();
              bool isInsideBone = false;
              double increase = -std::numeric_limits<double>::max() / points.size();
              unsigned numInBone = 0;
              for( int i = 0; i < points.size(); ++i) {
                  //statismo::VectorType normal = m_closestPoint->normalAtPoint(sample, m_tra)
                  //PointType closestPtOnSample =  m_closestPoint->findClosestPoint(sample, m_targetPoints[i]).first;
                  //double d = m_eval->evalSample(closestPtOnSample-m_targetPoints[i]);


                  if (m_closestPoint->huAtPoint(sample, i) > 1150 || m_closestPoint->huAtPoint(sample, i) < 400 )   {
//                      isInsideBone = true;
//                      break;

                      distance += increase;
                      numInBone +=  1;
                  }
                  else {
                      distance += 0;
                  }
              }

//              if (isInsideBone) {
//                  distance = -std::numeric_limits<double>::infinity();
//
//              } else {
//                  distance = 0;
//              }
//              std::cout << "number in bone is " << numInBone << std::endl;
              return distance;

          }
      private:
          const vector< PointType > m_targetPoints;
          const ActiveShapeModelType* m_asmodel;
          PositionEvaluator* m_eval;
          const ClosestPointType* m_closestPoint;

      };




      class Gaussian3DPositionDifferenceEvaluator : public PositionEvaluator
      {
        private:
          const RepresenterType* rep;
        public:
          Gaussian3DPositionDifferenceEvaluator( const RepresenterType* representer, double sigma = 1.0) : rep(representer){
            m_sigma = sigma;
            m_dNormalizer = -1.5*log( 2*M_PI ) - 3.0 * log( sigma );
          }


          double evalSample( PointType const& diff )
          {
            double dVal = rep->PointToVector(diff).squaredNorm(); // TODO check if PointType has this function norm2
            return -dVal/(2.0*m_sigma*m_sigma) + m_dNormalizer;
          }

          void evalGradient( PointType& grad, PointType const& diff )
          {
            //grad = diff*(-1) / (sigma * sigma);
          }

          virtual bool isSeparable() const
          {
            return true;
          }

        private:
          double m_sigma;
          double m_dNormalizer;
      };


      template <class T>
      class PointEvaluator : public DistributionEvaluator< MHFittingParameters > {
          typedef MeshOperations<typename Representer<T>::DatasetPointerType, typename Representer<T>::PointType> ClosestPointType;

        public:
          PointEvaluator( const Representer<T>* representer, const ClosestPointType* closestPoint, const vector< PointType >& targetPoints, ActiveShapeModelType* asmodel, PositionEvaluator* evaluator) :
            m_targetPoints(targetPoints),
            m_asmodel(asmodel),
            m_eval(evaluator),
            m_closestPoint(closestPoint)

          {}

          // DistributionEvaluatorInterface interface
        public:
          virtual double evalSample(const MHFittingParameters& currentSample) {
            double distance = 0.0;

            // TODO to draw full sample only for landmarks is inefficient

              typename RepresenterType::DatasetPointerType sample = m_closestPoint->transformMesh(currentSample);

            for( int i = 0; i < m_targetPoints.size(); ++i) {

              PointType closestPtOnSample =  m_closestPoint->findClosestPoint(sample, m_targetPoints[i]).first;
              double d = m_eval->evalSample(closestPtOnSample-m_targetPoints[i]);

                distance += d;
            }
            return distance;

          }
        private:
          const vector< PointType > m_targetPoints;
          const ActiveShapeModelType* m_asmodel;
          PositionEvaluator* m_eval;
          const ClosestPointType* m_closestPoint;

      };



      struct ASMLikelihoodForChunk {


          ASMLikelihoodForChunk(double _aggregatedLikelihood): aggregatedLikelihood(_aggregatedLikelihood) {  }

          double aggregatedLikelihood;

          // emulate move semantics, as boost::async seems to depend on it.
          ASMLikelihoodForChunk& operator=(BOOST_COPY_ASSIGN_REF(ASMLikelihoodForChunk) rhs) { // Copy assignment

              if (&rhs != this) {
                  copyMembers(rhs);
              }
              return *this;
          }

          ASMLikelihoodForChunk(BOOST_RV_REF(ASMLikelihoodForChunk) that) { //Move constructor
              copyMembers(that);
          }
          ASMLikelihoodForChunk& operator=(BOOST_RV_REF(ASMLikelihoodForChunk) rhs) { //Move assignment
              if (&rhs != this) {
                  copyMembers(rhs);
              }
              return *this;
          }
      private:
      BOOST_COPYABLE_AND_MOVABLE(ASMLikelihoodForChunk)
          void copyMembers(const ASMLikelihoodForChunk& that) {
              aggregatedLikelihood = that.aggregatedLikelihood;
          }
      };



/*
      template <class T>
      class ASMEvaluator : public DistributionEvaluator< MHFittingParameters > {
        typedef  ASMPreprocessedImage<typename RepresenterType::DatasetType> PreprocessedImageType;

      public:
          ASMEvaluator(ActiveShapeModelType* asmodel, PreprocessedImageType* image, PointSamplerType* sampler) :
                  m_asmodel(asmodel),
                  m_image(image),
                  m_sampler(sampler)
          {}

          // DistributionEvaluatorInterface interface
      public:
          virtual double evalSample(const MHFittingParameters& currentSample) {


              typename RepresenterType::DatasetPointerType sample = m_closestPoint->TransformMesh(m_asmodel->GetStatisticalModel(), currentSample);


              unsigned numProfilePointsUsed = 500;
              unsigned step = m_asmodel->GetProfiles().size() / numProfilePointsUsed;
//              FeatureExtractorType* fe = m_asmodel->GetFeatureExtractor()->CloneForTarget(m_asmodel,currentSample.GetCoefficients(),currentSample.GetRigidTransform());

              unsigned numChunks =  boost::thread::hardware_concurrency() + 1;
              std::vector<boost::future<ASMLikelihoodForChunk>* > futvec;


              for (unsigned i = 0; i < numChunks; ++i) {
                  unsigned profileIdStart = m_asmodel->GetProfiles().size() / numChunks  * i;
                  unsigned profileIdEnd = m_asmodel->GetProfiles().size() / numChunks * (i + 1);

                  boost::future<ASMLikelihoodForChunk> *fut = new boost::future<ASMLikelihoodForChunk>(
                          boost::async(boost::launch::async, boost::bind(&ASMEvaluator<T>::evalSampleForProfiles,
                                                                         this, profileIdStart, profileIdEnd, step,// fe,
                                                                         currentModelInstance, currentSample)));
                  futvec.push_back(fut);
              }


              double loglikelihood = 0.0;
              for (unsigned i = 0; i < futvec.size(); i++) {

                  ASMLikelihoodForChunk likelihoodForChunk = futvec[i]->get();
                  loglikelihood += likelihoodForChunk.aggregatedLikelihood;
                  delete futvec[i];
              }


              return loglikelihood;

          }


          ASMLikelihoodForChunk evalSampleForProfiles(unsigned profileIdStart, unsigned profileIdEnd, unsigned step, typename RepresenterType::DatasetPointerType currentModelInstance, const MHFittingParameters& currentSample) {

              FeatureExtractorType* fe = m_asmodel->GetFeatureExtractor()->CloneForTarget(m_asmodel,currentSample.GetCoefficients(),currentSample.GetRigidTransform());

              double loglikelihood = 0.0;

              for (unsigned i = profileIdStart; i < profileIdEnd; i += step) {

                  ASMProfile profile = m_asmodel->GetProfiles()[i];
                  long ptId = profile.GetPointId();

                  statismo::VectorType features;


                  bool ok = fe->ExtractFeatures(features, m_image, currentModelInstance->GetPoint(ptId));

                  if (ok) {
                      loglikelihood += profile.GetDistribution().logpdf(features);
                  } else {
                      std::cout << "feasutre not ok " << std::endl;
                  }
              }
              // evaluate profile points at ...
              fe->Delete();
//              ASMLikelihoodForChunk asmLikelihoodForChunk(loglikelihood);
              return ASMLikelihoodForChunk(loglikelihood);
          }

      private:
          const ActiveShapeModelType* m_asmodel;
          PreprocessedImageType* m_image;

          PointSamplerType* m_sampler;
          RandomGenerator* m_rGen;
      };
*/

      // TODO: Full image evaluator, i.e. ASM Evaluator is missing


    public:
      // "Script"
      /*============================================================================*/
      template <class T>
      class BasicSampling {
        public:
          static MarkovChain<MHFittingParameters >* buildChain(
                  const Representer<T>* representer,
                  const MeshOperations<typename Representer<T>::DatasetPointerType, typename Representer<T>::PointType>* closestPoint,
                  MHFittingConfiguration config,
                  ASMPreprocessedImage<typename Representer<T>::DatasetType>* targetImage,
              vector<PointType> targetPoints,
              ActiveShapeModelType* asmodel,
                  PointSamplerType* asmPointSampler,
              RigidTransformPointerType transform,
              statismo::VectorType coeffs) {

            // basics
            RandomGenerator* rGen = new RandomGenerator(42);
            MHFittingParameters init = MHFittingParameters(coeffs,transform);


            GaussianModelUpdate* poseAndModelProposalRough = new GaussianModelUpdate(0.2, 0.2, rGen);
            GaussianModelUpdate* poseAndModelProposalFine = new GaussianModelUpdate(0.05, 0.05, rGen);
              GaussianModelUpdate* poseAndModelProposalFinest = new GaussianModelUpdate(0.01, 0.01, rGen);
            //ASMModelUpdate* asmProposal = new ASMModelUpdate(config.GetAsmFittingconfiguration(), asmodel, targetImage, asmPointSampler);


            vector< typename RandomProposal< MHFittingParameters >::GeneratorPair> gaussMixtureProposalVector(3);
            gaussMixtureProposalVector[0] = pair<ProposalGenerator<MHFittingParameters >*,double>(poseAndModelProposalRough,0.02);
            gaussMixtureProposalVector[1] = pair<ProposalGenerator<MHFittingParameters >*,double>(poseAndModelProposalFine,0.2);
              gaussMixtureProposalVector[2] = pair<ProposalGenerator<MHFittingParameters >*,double>(poseAndModelProposalFinest,0.6);

              RandomProposal<MHFittingParameters >* gaussMixtureProposal = new RandomProposal<MHFittingParameters >(gaussMixtureProposalVector,rGen);

            Gaussian3DPositionDifferenceEvaluator* diffEval = new Gaussian3DPositionDifferenceEvaluator(asmodel->GetRepresenter(), 0.4);
            PointEvaluator<T>* pointEval = new PointEvaluator<T>(representer, closestPoint, targetPoints,asmodel,diffEval);
            ModelPriorEvaluator* modelPriorEvaluator = new ModelPriorEvaluator();
//

              InLungOrBoneEvaluator<T>* huEvaluator = new InLungOrBoneEvaluator<T>(representer, closestPoint, asmodel);

              std::vector<DistributionEvaluator<MHFittingParameters >*> huEvaluatorList;
              huEvaluatorList.push_back(huEvaluator);
              huEvaluatorList.push_back(modelPriorEvaluator);


              QuietLogger <MHFittingParameters>* ql = new QuietLogger<MHFittingParameters>();
              MarkovChain<MHFittingParameters >* huChain =  new MetropolisHastings<MHFittingParameters >(gaussMixtureProposal, new ProductEvaluator<MHFittingParameters>(huEvaluatorList), ql, init, rGen );
              MarkovChainProposal<MHFittingParameters>* huChainProposal = new MarkovChainProposal<MHFittingParameters>(huChain, 10);




              std::vector<DistributionEvaluator<MHFittingParameters >*> lmAndHuEvaluatorList;
              lmAndHuEvaluatorList.push_back(pointEval);
              lmAndHuEvaluatorList.push_back(huEvaluator);
              lmAndHuEvaluatorList.push_back(modelPriorEvaluator);

              MarkovChain<MHFittingParameters >* lmChain =  new MetropolisHastings<MHFittingParameters >(huChainProposal, new ProductEvaluator<MHFittingParameters>(lmAndHuEvaluatorList), ql, init, rGen );
              MarkovChainProposal<MHFittingParameters>* huAndLMProposal = new MarkovChainProposal<MHFittingParameters>(lmChain, 10);


//            Gaussian3DPositionDifferenceEvaluator* wideDiffEval = new Gaussian3DPositionDifferenceEvaluator(asmodel->GetRepresenter(), 5.0);
//            PointEvaluator<T>* widePointEval = new PointEvaluator<T>(representer, closestPoint, targetPoints,asmodel,wideDiffEval);
            //ASMEvaluator<T>* asmEvaluator = new ASMEvaluator<T>(asmodel, targetImage, asmPointSampler);

              //std::vector<DistributionEvaluator<MHFittingParameters >*> completeEvaluatorList;
              //vector< typename RandomProposal< MHFittingParameters >::GeneratorPair> finalProposalVector;


              //completeEvaluatorList.push_back(pointEval);
              //completeEvaluatorList.push_back(asmEvaluator);
              //completeEvaluatorList.push_back(huEvaluator);
              //completeEvaluatorList.push_back(modelPriorEvaluator);

              //finalProposalVector.push_back(pair<ProposalGenerator<MHFittingParameters >*,double>(huAndLMProposal,0.99));
              //finalProposalVector.push_back(pair<ProposalGenerator<MHFittingParameters >*,double>(asmProposal,0.01));


              //RandomProposal<MHFittingParameters >* finalProposal = new RandomProposal<MHFittingParameters >(finalProposalVector,rGen);



            // markov chain
              //BestMatchLogger<MHFittingParameters>* bml = new BestMatchLogger<MHFittingParameters>();

              // WARNING!  we currently do not use the ASM proposal
              MarkovChain<MHFittingParameters >* chain = new MetropolisHastings<MHFittingParameters>(huAndLMProposal, new ProductEvaluator<MHFittingParameters >(lmAndHuEvaluatorList), ql, init, rGen);

            return chain;
          }

          static MarkovChain<MHFittingParameters >* buildLmAndHuChain(
            const Representer<T>* representer,
            const MeshOperations<typename Representer<T>::DatasetPointerType, typename Representer<T>::PointType>* closestPoint,
            vector<PointType> targetPoints,
            ActiveShapeModelType* asmodel,
            RigidTransformPointerType transform,
            statismo::VectorType coeffs) {

            // basics
            RandomGenerator* rGen = new RandomGenerator(42);
            MHFittingParameters init = MHFittingParameters(coeffs, transform);

            GaussianModelUpdate* poseAndShapeProposalRough = new GaussianModelUpdate(0.5, 0.5, rGen);
            GaussianModelUpdate* poseAndShapeProposalFine = new GaussianModelUpdate(0.2, 0.2, rGen);
            GaussianModelUpdate* poseAndShapeProposalFinest = new GaussianModelUpdate(0.05, 0.05, rGen);


            vector< typename RandomProposal< MHFittingParameters >::GeneratorPair> gaussMixtureProposalVector(3);
            gaussMixtureProposalVector[0] = pair<ProposalGenerator<MHFittingParameters >*, double>(poseAndShapeProposalRough, 0.02);
            gaussMixtureProposalVector[1] = pair<ProposalGenerator<MHFittingParameters >*, double>(poseAndShapeProposalFine, 0.2);
            gaussMixtureProposalVector[2] = pair<ProposalGenerator<MHFittingParameters >*, double>(poseAndShapeProposalFinest, 0.6);

            RandomProposal<MHFittingParameters >* gaussMixtureProposal = new RandomProposal<MHFittingParameters >(gaussMixtureProposalVector, rGen);

            Gaussian3DPositionDifferenceEvaluator* diffEval = new Gaussian3DPositionDifferenceEvaluator(asmodel->GetRepresenter(), 0.4);
            PointEvaluator<T>* pointEval = new PointEvaluator<T>(representer, closestPoint, targetPoints, asmodel, diffEval);
            ModelPriorEvaluator* modelPriorEvaluator = new ModelPriorEvaluator();

            InLungOrBoneEvaluator<T>* huEvaluator = new InLungOrBoneEvaluator<T>(representer, closestPoint, asmodel);

            std::vector<DistributionEvaluator<MHFittingParameters >*> huEvaluatorList;
            huEvaluatorList.push_back(huEvaluator);
            huEvaluatorList.push_back(modelPriorEvaluator);
            
            QuietLogger <MHFittingParameters>* ql = new QuietLogger<MHFittingParameters>();
            MarkovChain<MHFittingParameters >* huChain = new MetropolisHastings<MHFittingParameters >(gaussMixtureProposal, new ProductEvaluator<MHFittingParameters>(huEvaluatorList), ql, init, rGen);
            MarkovChainProposal<MHFittingParameters>* huChainProposal = new MarkovChainProposal<MHFittingParameters>(huChain, 10);
            


            std::vector<DistributionEvaluator<MHFittingParameters >*> lmAndHuEvaluatorList;
            lmAndHuEvaluatorList.push_back(pointEval);
            lmAndHuEvaluatorList.push_back(huEvaluator);
            lmAndHuEvaluatorList.push_back(modelPriorEvaluator);

            MarkovChain<MHFittingParameters >* lmAndHuChain = new MetropolisHastings<MHFittingParameters >(huChainProposal, new ProductEvaluator<MHFittingParameters>(lmAndHuEvaluatorList), ql, init, rGen);
            MarkovChainProposal<MHFittingParameters>* lmAndHuProposal = new MarkovChainProposal<MHFittingParameters>(lmAndHuChain, 10);
            
            MarkovChain<MHFittingParameters >* chain = new MetropolisHastings<MHFittingParameters>(lmAndHuProposal, new ProductEvaluator<MHFittingParameters >(lmAndHuEvaluatorList), ql, init, rGen);

            return chain;
          }

          static MarkovChain<MHFittingParameters >* buildInitialPoseChain(
            const Representer<T>* representer,
            const MeshOperations<typename Representer<T>::DatasetPointerType, typename Representer<T>::PointType>* closestPoint,
            vector<PointType> targetPoints,
            ActiveShapeModelType* asmodel,
            MHFittingParameters& initialParameters) {

            // basics
            RandomGenerator* rGen = new RandomGenerator(42);


            GaussianModelUpdate* poseProposalRough = new GaussianModelUpdate(0, 1, rGen);
            GaussianModelUpdate* poseProposalFine = new GaussianModelUpdate(0, 0.3, rGen);
            GaussianModelUpdate* poseProposalFinest = new GaussianModelUpdate(0, 0.1, rGen);


            vector< typename RandomProposal< MHFittingParameters >::GeneratorPair> gaussMixtureProposalVector(3);
            gaussMixtureProposalVector[0] = pair<ProposalGenerator<MHFittingParameters >*, double>(poseProposalRough, 10);
            gaussMixtureProposalVector[1] = pair<ProposalGenerator<MHFittingParameters >*, double>(poseProposalFine, 1);
            gaussMixtureProposalVector[2] = pair<ProposalGenerator<MHFittingParameters >*, double>(poseProposalFinest, 0.1);

            RandomProposal<MHFittingParameters >* gaussMixtureProposal = new RandomProposal<MHFittingParameters >(gaussMixtureProposalVector, rGen);

            Gaussian3DPositionDifferenceEvaluator* diffEval = new Gaussian3DPositionDifferenceEvaluator(asmodel->GetRepresenter(), 0.4);
            PointEvaluator<T>* pointEval = new PointEvaluator<T>(representer, closestPoint, targetPoints, asmodel, diffEval);
            
            QuietLogger <MHFittingParameters>* ql = new QuietLogger<MHFittingParameters>();
            MarkovChain<MHFittingParameters >* lmChain = new MetropolisHastings<MHFittingParameters >(gaussMixtureProposal, pointEval, ql, initialParameters, rGen);

            return lmChain;
          }


          static MarkovChain<MHFittingParameters >* buildInitialModelChain(
            const Representer<T>* representer,
            const MeshOperations<typename Representer<T>::DatasetPointerType, typename Representer<T>::PointType>* closestPoint,
            vector<PointType> targetPoints,
            ActiveShapeModelType* asmodel,
            MHFittingParameters& initialParameters) {

            // basics
            RandomGenerator* rGen = new RandomGenerator(42);


            GaussianModelUpdate* poseAndShapeProposalRough = new GaussianModelUpdate(1, 1, rGen);
            GaussianModelUpdate* poseAndShapeProposalFine = new GaussianModelUpdate(0.3, 0.3, rGen);
            GaussianModelUpdate* poseAndShapeProposalFinest = new GaussianModelUpdate(0.1, 0.1, rGen);


            vector< typename RandomProposal< MHFittingParameters >::GeneratorPair> gaussMixtureProposalVector(3);
            gaussMixtureProposalVector[0] = pair<ProposalGenerator<MHFittingParameters >*, double>(poseAndShapeProposalRough, 0.02);
            gaussMixtureProposalVector[1] = pair<ProposalGenerator<MHFittingParameters >*, double>(poseAndShapeProposalFine, 0.2);
            gaussMixtureProposalVector[2] = pair<ProposalGenerator<MHFittingParameters >*, double>(poseAndShapeProposalFinest, 0.6);

            RandomProposal<MHFittingParameters >* gaussMixtureProposal = new RandomProposal<MHFittingParameters >(gaussMixtureProposalVector, rGen);

            Gaussian3DPositionDifferenceEvaluator* diffEval = new Gaussian3DPositionDifferenceEvaluator(asmodel->GetRepresenter(), 0.4);
            PointEvaluator<T>* pointEval = new PointEvaluator<T>(representer, closestPoint, targetPoints, asmodel, diffEval);

            ModelPriorEvaluator* modelPriorEvaluator = new ModelPriorEvaluator();

            std::vector<DistributionEvaluator<MHFittingParameters >*> lmAndPriorEvaluatorList;
            lmAndPriorEvaluatorList.push_back(pointEval);
            lmAndPriorEvaluatorList.push_back(modelPriorEvaluator);

            QuietLogger <MHFittingParameters>* ql = new QuietLogger<MHFittingParameters>();
            MarkovChain<MHFittingParameters >* lmChain = new MetropolisHastings<MHFittingParameters >(gaussMixtureProposal, new ProductEvaluator<MHFittingParameters >(lmAndPriorEvaluatorList), ql, initialParameters, rGen);

            return lmChain;
          }





      };

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

      MHFittingStep(MarkovChain<MHFittingParameters >* chain)
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


      static MHFittingStep *Create(MarkovChain<MHFittingParameters >* chain) {
        return new MHFittingStep( chain);
      }

      ~MHFittingStep() {
      }

      MHFittingParameters Perform() const {


        //          ASMFittingStepType* asmFittingStep = ASMFittingStepType::Create(m_configuration.GetAsmFittingconfiguration(), m_model, m_sourceCoefficients, m_sourceTransform, m_target, m_sampler);
        //          ASMFittingResult<RigidTransformPointerType> result = asmFittingStep->Perform();


        // runs a markov chain and returns only the accepted proposals
        MHFittingParameters params;//(m_sourceTransform,m_sourceCoefficients);
        m_chain->next(params);
        MHFittingParameters mhResult(params.GetCoefficients(), params.GetRigidTransformParameters());
        return mhResult;
      }

    private:
      MarkovChain<MHFittingParameters >* m_chain;



      statismo::VectorType fromVnlVector(const vnl_vector<float>& v) {
          return Eigen::Map<const statismo::VectorType>(v.data_block(), v.size());

      }


  };
}

#endif //STATISMO_ASMFITTING_H


