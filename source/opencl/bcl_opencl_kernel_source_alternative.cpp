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
#include "opencl/bcl_opencl_kernel_source_alternative.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from two alternatives
    //! @param ALTERNATIVE_SOURCE_1
    //! @param ALTERNATIVE_SOURCE_2
    KernelSourceAlternative::KernelSourceAlternative
    (
      const KernelSourceInterface &ALTERNATIVE_SOURCE_1,
      const KernelSourceInterface &ALTERNATIVE_SOURCE_2
    ) :
      m_Alternatives()
    {
      m_Alternatives.PushBack( util::ShPtr< KernelSourceInterface>( ALTERNATIVE_SOURCE_1.Clone()));
      m_Alternatives.PushBack( util::ShPtr< KernelSourceInterface>( ALTERNATIVE_SOURCE_2.Clone()));
    }

    //! @brief Clone function
    //! @return pointer to new KernelSourceAlternative
    KernelSourceAlternative *KernelSourceAlternative::Clone() const
    {
      return new KernelSourceAlternative( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &KernelSourceAlternative::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief identifier for kernel source, so that the progam build can be cached based on it
    //! @return identifier like the filename or the function name
    const std::string &KernelSourceAlternative::GetIdentifier() const
    {
      static const std::string s_identifier;
      return m_Alternatives.IsEmpty() ? s_identifier : m_Alternatives.FirstElement()->GetIdentifier();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get the source for compilation
    //! @param PRECISION precision of the kernel
    //! @param EXTENSIONS extensions of the devices, this kernel is compiled for
    //! @return source with precision set correctly
    std::string KernelSourceAlternative::GetSource( const util::CPPDataTypes::Types PRECISION, const storage::Set< Extension> &EXTENSIONS) const
    {
      // iterate over all alternatives
      for
      (
        util::ShPtrVector< KernelSourceInterface>::const_iterator
          itr( m_Alternatives.Begin()), itr_end( m_Alternatives.End());
        itr != itr_end;
        ++itr
      )
      {
        const std::string source( ( *itr)->GetSource( PRECISION, EXTENSIONS));
        if( !source.empty())
        {
          return source;
        }
      }

      // no alternative gave source
      return std::string();
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &KernelSourceAlternative::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Alternatives, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &KernelSourceAlternative::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Alternatives, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace opencl
} // namespace bcl
