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
#include "biol/bcl_biol_mutation.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "biol/bcl_biol_atom.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_string_functions.h"
#include "util/bcl_util_string_numeric_conversion.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new Mutation
    Mutation *Mutation::Clone() const
    {
      return new Mutation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Mutation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief write the mutation as a string in standard mutation format
    std::string Mutation::ToString() const
    {
      return std::string( size_t( 1), m_NativeType->GetOneLetterCode())
             + util::Format()( m_ResidueNumber)
             + std::string( size_t( 1), m_MutantType->GetOneLetterCode());
    }

    //! @brief write the mutation as a string in standard mutation format
    Mutation Mutation::FromString( const std::string &STR)
    {
      BCL_Assert( STR.length() >= size_t( 3), "need at least 3 characters for a valid mutation!");
      size_t end( STR.find_first_not_of( "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"));
      if( end == std::string::npos)
      {
        end = STR.size();
      }
      return
        Mutation
        (
          util::ConvertStringToNumericalValue< int>( STR.substr( 1, end - 2)),
          GetAATypes().AATypeFromOneLetterCode( STR[0]),
          GetAATypes().AATypeFromOneLetterCode( STR[end -1]),
          end + size_t( 1) < STR.size() ? STR.substr( end + size_t( 1)) : ""
        );
    }

    //! @param WITH_DATA whether to include any data members, else, only include initialization members
    util::ObjectDataLabel Mutation::GetLabel( const bool &WITH_DATA) const
    {
      return util::ObjectDataLabel( "", ToString());
    }

    //! @brief connect the mutation with a particular residue
    void Mutation::SetAA( const AABase &BASE) const
    {
      for( auto itr( m_AAs.Begin()), itr_end( m_AAs.End()); itr != itr_end; ++itr)
      {
        if( ( *itr)->GetSeqID() == BASE.GetSeqID() && ( *itr)->GetChainID() == BASE.GetChainID())
        {
          *itr = util::ToSiPtr( BASE);
          return;
        }
      }
      m_AAs.PushBack( util::ToSiPtr( BASE));
      if( m_ChainIDs.find( BASE.GetChainID()) == std::string::npos)
      {
        m_ChainIDs += BASE.GetChainID();
      }
      BCL_MessageVrb( "New aa: " + BASE.GetIdentification() + " " + util::Format()( BASE.GetCA().GetCoordinates()));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief set the value of the corresponding member based on the label
    //! @param LABEL label that is used to set the string
    //! @param ERROR_STREAM stream to write errors to
    //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
    bool Mutation::TryRead( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      if( LABEL.GetValue().size() < size_t( 3))
      {
        this->WriteHelp( ERROR_STREAM);
        return false;
      }
      m_NativeType = GetAATypes().AATypeFromOneLetterCode( LABEL.GetValue()[ 0]);
      if( !m_NativeType.IsDefined())
      {
        ERROR_STREAM << LABEL.GetValue()[ 0] << " in mutation " << LABEL.ToString() << " is not a valid amino acid type";
        return false;
      }
      size_t end( LABEL.GetValue().find_first_not_of( "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"));
      if( end == std::string::npos)
      {
        end = LABEL.GetValue().size();
      }
      m_MutantType = GetAATypes().AATypeFromOneLetterCode( LABEL.GetValue()[ end - size_t( 1)]);
      if( !m_MutantType.IsDefined())
      {
        ERROR_STREAM << LABEL.GetValue()[ LABEL.GetValue().size() - size_t( 1)] << " in mutation "
                     << LABEL.ToString() << " is not a valid amino acid type";
        return false;
      }
      if( !util::TryConvertFromString( m_ResidueNumber, LABEL.GetValue().substr( 1, end - 2), ERROR_STREAM))
      {
        ERROR_STREAM << " in mutation " << LABEL.ToString();
        return false;
      }
      m_ChainIDs = end + size_t( 1) < LABEL.GetValue().size() ? LABEL.GetValue().substr( end + size_t( 1)) : "";
      return true;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Mutation::Read( std::istream &ISTREAM)
    {
      BCL_Assert( this->TryRead( util::ObjectDataLabel( ISTREAM), util::GetLogger()), "Could not read mutation");
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &Mutation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM << ToString();
    }

    //! @brief writes the help for the label
    //! @param OSTREAM the stream to which the help is written to
    //! @param INDENT the amount of indent
    //! @return the given stream to which the help was written to
    std::ostream &Mutation::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM << "Mutation in standard format: NativeAATypeResIDMutantAAType like V215M. May optionally indicate chain"
                     << " id(s) afterwards, e.g. V215M-ABCD";
    }

  } // namespace biol
} // namespace bcl
