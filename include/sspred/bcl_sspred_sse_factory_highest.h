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

#ifndef BCL_SSPRED_SSE_FACTORY_HIGHEST_H_
#define BCL_SSPRED_SSE_FACTORY_HIGHEST_H_

// include the namespace header
#include "bcl_sspred.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sspred_methods.h"
#include "assemble/bcl_assemble_sse_factory_interface.h"
#include "biol/bcl_biol_ss_types.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEFactoryHighest
    //! @brief an SSEFactory implementation that uses the highest of all predictions to assign secondary structure
    //! @details
    //!
    //! @see @link example_sspred_sse_factory_highest.cpp @endlink
    //! @author woetzen
    //! @date Jun 15, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEFactoryHighest :
      public assemble::SSEFactoryInterface
    {

    private:

    //////////
    // data //
    //////////

      //! secondary structure prediction method to use
      Method m_Method;

      //! threshold for pool ssprediction structure sstype assignment
      storage::Map< biol::SSType, double> m_MinStructureThreshold;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from sspred method
      //! @param SSMETHOD sspred method to use to generate pool of sses
      SSEFactoryHighest
      (
        const Method &SSMETHOD
      );

      //! @brief Clone function
      //! @return pointer to new SSEFactoryHighest
      SSEFactoryHighest *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that returns a set of SSEs for the given AASequence
      //! @brief SEQUENCE AASequence from which the SSEPool is going to be built
      //! @return SSEPool built from provided SEQUENCE
      assemble::SSEPool
      operator()( const biol::AASequence &SEQUENCE) const;

      //! @brief set the ss method used
      //! @param SS_METHOD sspred method
      void SetMethod( const Method &SS_METHOD)
      {
        m_Method = SS_METHOD;
      }

      //! @brief set the thresholds to use
      //! @brief SSTYPE_THRESHOLDS the thresholds to use for the sstypes desired
      void SetThresholds( const storage::Map< biol::SSType, double> &SSTYPE_THRESHOLDS)
      {
        m_MinStructureThreshold = SSTYPE_THRESHOLDS;
      }

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

    }; // class SSEFactoryHighest

  } // namespace sspred
} // namespace bcl

#endif // BCL_SSPRED_SSE_FACTORY_HIGHEST_H_ 
