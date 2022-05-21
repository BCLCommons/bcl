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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "sspred/bcl_sspred_sse_factory_highest.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_pool.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEFactoryHighest::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEFactoryHighest( GetMethods().e_Undefined))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from sspred method
    //! @param SSMETHOD sspred method to use to generate pool of sses
    SSEFactoryHighest::SSEFactoryHighest
    (
      const Method &SSMETHOD
    ) :
      m_Method( SSMETHOD)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSEFactoryHighest
    SSEFactoryHighest *SSEFactoryHighest::Clone() const
    {
      return new SSEFactoryHighest( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEFactoryHighest::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that returns a set of SSEs for the given AASequence
    //! @brief SEQUENCE AASequence from which the SSEPool is going to be built
    //! @return SSEPool built from provided SEQUENCE
    assemble::SSEPool
    SSEFactoryHighest::operator()( const biol::AASequence &SEQUENCE) const
    {
      // initialize a new pool
      assemble::SSEPool sse_pool;

      // initialize a variable to hold the last type of the SSE generated
      biol::SSType sstype;

      // make a copy of the sequence and reset
      biol::AASequence sequence( SEQUENCE);
      sequence.Reset();

      // iterate over sequence
      for( biol::AASequence::const_iterator itr( SEQUENCE.Begin()), itr_end( SEQUENCE.End()); itr != itr_end; ++itr)
      {
        // get one state of the current AA from its predictions
        biol::SSType current_sstype( ( *itr)->GetSSPrediction( m_Method)->GetOneStateSSPrediction());

        // find the corresponding threshold for the given type
        const storage::Map< biol::SSType, double>::const_iterator min_th_itr( m_MinStructureThreshold.Find( current_sstype));

        // if given type is not a coil and threshold was found and the value is larger smaller than threshold
        if
        (
          current_sstype != biol::GetSSTypes().COIL &&
          min_th_itr != m_MinStructureThreshold.End() &&
          ( *itr)->GetSSPrediction( m_Method)->GetThreeStatePrediction()( current_sstype) < min_th_itr->second
        )
        {
          // update sstype to be coil
          current_sstype = biol::GetSSTypes().COIL;
        }

        // insert the previous sequence as sse if the type of the current aa changes
        if( current_sstype != sstype)
        {
          if( sequence.GetSize() > 0)
          {
            util::ShPtr< assemble::SSE> sp_sse( new assemble::SSE( sequence, sstype));
            sse_pool.Insert( sp_sse);
            sequence.Reset();
          }

          // update sstype and reset
          sstype = current_sstype;
        }
        // insert amino acid into sequence
        sequence.PushBack( *itr);
      }

      // insert the remaining sequence
      if( sequence.GetSize() > 0)
      {
        util::ShPtr< assemble::SSE> sp_sse( new assemble::SSE( sequence, sstype));
        sse_pool.Insert( sp_sse);
      }

      // end
      return sse_pool;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEFactoryHighest::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Method               , ISTREAM);
      io::Serialize::Read( m_MinStructureThreshold, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEFactoryHighest::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Method               , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinStructureThreshold, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace sspred
} // namespace bcl
