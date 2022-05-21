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

#ifndef BCL_ASSEMBLE_SSE_FACTORY_CONFORMATION_H_
#define BCL_ASSEMBLE_SSE_FACTORY_CONFORMATION_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_factory_interface.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEFactoryConformation
    //! @brief This is an SSEFactory class that uses backbone conformation to identify SSEs for a given AASequence
    //! @details This class is derived from SSEFactoryInterface and it uses real backbone conformations, specifically
    //! phi-psi angle to determine SSEs for a given AASequence. This class is important for getting a ProteinModel
    //! with SSEs when a given PDB files lacks the SSE definitions
    //!
    //! @see @link example_assemble_sse_factory_conformation.cpp @endlink
    //! @author woetzen, karakam
    //! @date 20.03.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEFactoryConformation :
      public SSEFactoryInterface
    {

    private:

    //////////
    // data //
    //////////

      //! minimal sequence length to start angle calculation
      static const size_t s_MinimalSequenceLength = 3;

      //! minimum sizes for each SSType
      storage::VectorND< 3, size_t> m_MinSSELengths;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEFactoryConformation();

      //! @brief constructor from provided minimal sse lengths
      //! @param MIN_SSE_LENGTHS minimum size of SSEs to be considered
      SSEFactoryConformation( const storage::VectorND< 3, size_t> &MIN_SSE_LENGTHS);

      //! @brief virtual copy constructor
      virtual SSEFactoryConformation *Clone() const;

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
      }

      //! @brief set the thresholds to use
      //! @brief SSTYPE_THRESHOLDS the thresholds to use for the sstypes desired
      void SetThresholds( const storage::Map< biol::SSType, double> &SSTYPE_THRESHOLDS)
      {
      }

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
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class SSEFactoryConformation

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_SSE_FACTORY_CONFORMATION_H_
