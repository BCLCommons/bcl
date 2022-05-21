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
#include "fold/bcl_fold_mutate_sse_type.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateSSEType::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateSSEType( storage::Set< biol::SSType>()))
    );

    //! @brief the default scheme
    //! @return reference to default scheme as string
    const std::string &MutateSSEType::GetDefaultScheme()
    {
      static const std::string s_default_scheme( "sse_type");
      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! @param POSSIBLE_TYPES set of types that sses can be set to
    //! @param SCHEME the scheme of the mutate
    MutateSSEType::MutateSSEType
    (
      const storage::Set< biol::SSType> &POSSIBLE_TYPES,
      const std::string &SCHEME
    ) :
      m_PossibleTypes( POSSIBLE_TYPES),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateSSEType
    MutateSSEType *MutateSSEType::Clone() const
    {
      return new MutateSSEType( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateSSEType::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief change the sse type and the conformation accordingly
    //! @param ELEMENT the sse to change
    //! @return MutateResult that results from mutating to the SSE
    math::MutateResult< assemble::SSE> MutateSSEType::operator()( const assemble::SSE &ELEMENT) const
    {
      // reduce the set of types to be considered
      storage::Set< biol::SSType> remaining_types( m_PossibleTypes);
      remaining_types.Erase( ELEMENT.GetType());

      // if there are no types left to set it to
      if( remaining_types.IsEmpty())
      {
        return math::MutateResult< assemble::SSE>( util::ShPtr< assemble::SSE>(), *this);
      }

      // get a random new type
      const storage::Set< biol::SSType>::const_iterator itr
      (
        random::GetGlobalRandom().Iterator( remaining_types.Begin(), remaining_types.End(), remaining_types.GetSize())
      );

      // create a new sse from the sequence and the new sstype
      const util::ShPtr< assemble::SSE> sp_sse( new assemble::SSE( ELEMENT, *itr));

      BCL_MessageDbg( "changed sstype: " + ELEMENT.GetIdentification() + '\t' + sp_sse->GetIdentification());

      // return the result
      return math::MutateResult< assemble::SSE>( sp_sse, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateSSEType::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_PossibleTypes, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateSSEType::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_PossibleTypes, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme,        OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
