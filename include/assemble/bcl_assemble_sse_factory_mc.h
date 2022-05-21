// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_ASSEMBLE_SSE_FACTORY_MC_H_
#define BCL_ASSEMBLE_SSE_FACTORY_MC_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "sspred/bcl_sspred.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_factory_interface.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "math/bcl_math_mutate_interface.h"
#include "util/bcl_util_enum.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEFactoryMC
    //! @brief an SSEFactory implementation that uses the SSEConfidence score to derive a pool of sselements
    //! @details using a monte carlo minimization scheme, the sse confidence over the pool is optimized, while the
    //!          aa sequence is split into secondary structure elements and amino acids are moved from one to another
    //!          sse. A single secondary structure prediction method is used.
    //!
    //! @see @link example_assemble_sse_factory_mc.cpp @endlink
    //! @author woetzen
    //! @date Jun 15, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEFactoryMC :
      public SSEFactoryInterface
    {

    private:

    //////////
    // data //
    //////////

      //! secondary structure prediction method to use
      sspred::Method m_Method;

      //! confidence threshold; z-score threshold above which the confidence score is negative
      double m_ConfidenceThreshold;

      //! scoring function used in minimization
      util::ShPtr< math::BinaryFunctionInterface< SSEPool, biol::Membrane, double> > m_ScoringFunction;

      //! mutate used inerate minimization
      util::ShPtr< math::MutateInterface< SSEPool> > m_Mutate;

      //! max number of iterations
      size_t m_MaxNumberIterations;

      //! number of optimizations, of which the best one is returned
      size_t m_NumberOptimizations;

      //! Defined constants / default values
      enum
      {
        s_MaxNumberIterations = 1000, //!< Default maximum number of iterations
        s_MaxNumberUnimproved =  250, //!< Default max number of unimproved steps
        s_NumberOptimizations =    1  //!< Default number of optimizations, of which the best one is returned
      };

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from sspred method
      //! @param SSMETHOD sspred method to use to generate pool of sses
      //! @param CONFIDENCE_THRESHOLD the threshold in units of z-score, above which the confidence score becomes negative
      SSEFactoryMC
      (
        const sspred::Method &SSMETHOD,
        const double CONFIDENCE_THRESHOLD
      );

      //! @brief Clone function
      //! @return pointer to new SSEFactoryMC
      SSEFactoryMC *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief set the ss method used
      //! @param SS_METHOD sspred method
      void SetMethod( const sspred::Method &SS_METHOD)
      {
        m_Method = SS_METHOD;
      }

      //! @brief set the thresholds to use
      //! @brief SSTYPE_THRESHOLDS the thresholds to use for the sstypes desired
      void SetThresholds( const storage::Map< biol::SSType, double> &SSTYPE_THRESHOLDS);

      //! @brief set the scoring function
      //! @param SP_SCORING ShPtr to scoring function
      void SetScoringFunction( const util::ShPtr< math::BinaryFunctionInterface< SSEPool, biol::Membrane, double> > &SP_SCORING);

      //! @brief access to the scoring function
      //! @return the scoring function used in minimization
      const util::ShPtr< math::BinaryFunctionInterface< SSEPool, biol::Membrane, double> > &GetScoringFunction() const
      {
        return m_ScoringFunction;
      }

      //! @brief set the mutate object
      //! @param SP_MUTATE ShPtr to mutate
      void SetMutate( const util::ShPtr< math::MutateInterface< SSEPool> > &SP_MUTATE);

      //! @brief set the number of optimizations to perform, of which the best is the pool beeing returned by the operator
      //! @param NUMBER_OPTIMZATIONS
      void SetNumberOptimizations( const size_t NUMBER_OPTIMZATIONS)
      {
        m_NumberOptimizations = NUMBER_OPTIMZATIONS;
      }

      //! @brief set the max number of iterations per minimization
      //! @param MAX_NUMBER_ITERATIONS maximal number of MC steps
      void SetMaxNumberIterations( const size_t MAX_NUMBER_ITERATIONS)
      {
        m_MaxNumberIterations = MAX_NUMBER_ITERATIONS;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief generate the default scoring function
      //! @return function to score the sse pool as objective function
      util::ShPtr< math::BinaryFunctionInterface< SSEPool, biol::Membrane, double> > DefaultScoringFunction() const;

      //! @brief generate the default set of mutates
      //! @return mutate to change the compositon of the sse pool
      util::ShPtr< math::MutateInterface< SSEPool> > DefaultMutate() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that returns a set of SSEs for the given AASequence
      //! @brief SEQUENCE AASequence from which the SSEPool is going to be built
      //! @return SSEPool built from provided SEQUENCE
      SSEPool
      operator()( const biol::AASequence &SEQUENCE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class SSEFactoryMC

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_FACTORY_MC_H_
