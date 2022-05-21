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
#include "nmr/bcl_nmr_signal_1d.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {

    //! @brief PatternType as string
    //! @param PATTERN_TYPE the PatternType
    //! @return the string for the PatternType
    const std::string &Signal1D::GetPatternTypeDescriptor( const PatternType &PATTERN_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "UnknownPattern",
        "Singlet",
        "Doublet",
        "Triplet",
        "Quadruplet",
        "Multiplet",
        GetStaticClassName< PatternType>()
      };

      return s_descriptors[ PATTERN_TYPE];
    }

    //! @brief PatternType as single character
    //! @details used in mdl Spectrum lines to describe the pattern
    //! @param PATTERN_TYPE the PatternType
    //! @return the char for the PatternType
    const char &Signal1D::GetPatternTypeSingleLetterCode( const PatternType &PATTERN_TYPE)
    {
      static const std::string s_single_letter_code( "USDTQMX");

      return s_single_letter_code[ PATTERN_TYPE];
    }

    //! @brief PatternType from single letter code
    //! @details convert char into a PatternType
    //! @param SINGLE_LETTER_CODE the single letter code as it may appear in an mdl Spectrum line
    //! @return PatternType for single letter code, e_UnknownPattern if not found
    Signal1D::PatternType Signal1D::PatternTypeFromSingleLetterCode( const char SINGLE_LETTER_CODE)
    {
      // range to search
      const char *beg( &GetPatternTypeSingleLetterCode( e_UnknownPattern));
      const char *end( beg + s_NumberPatternTypes);

      // find the single letter code
      const char *pattern( std::find( beg, end, SINGLE_LETTER_CODE));

      // if it was not found, return unknown
      if( pattern == end)
      {
        return e_UnknownPattern;
      }

      // return the PatternType
      return PatternType( pattern - beg);
    }

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Signal1D::s_Instance
    (
      GetObjectInstances().AddInstance( new Signal1D())
    );

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Signal1D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read Signal1D from std::istream
    std::istream &Signal1D::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ChemicalShift, ISTREAM);
      io::Serialize::Read( m_Integral,      ISTREAM);
      io::Serialize::Read( m_Pattern,       ISTREAM);
      io::Serialize::Read( m_Atom,          ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write Signal1D into std::ostream
    std::ostream &Signal1D::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ChemicalShift, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Integral,      OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Pattern,       OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Atom,          OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Signal from single entry in mdl
    //! @param SIGNAL_STRING string of the form shift;coupling_constcoupling;atom_index
    //! @param ATOMS the atoms, the atom index is relative to
    //! @return Signal1D for that string
    util::ShPtr< Signal1D> Signal1D::SignalFromString
    (
      const std::string &SIGNAL_STRING,
      iterate::Generic< const chemistry::AtomConformationalInterface> ATOMS
    )
    {
      // split the signal string
      const storage::Vector< std::string> signal_split( util::SplitString( SIGNAL_STRING, std::string( 1, s_InformationSeparator)));
      if( signal_split.GetSize() != 3)
      {
        BCL_MessageCrt
        (
          "each individual signal needs to contain: shift;coupling;atim_index, but it is: " + SIGNAL_STRING
        );
        return util::ShPtr< Signal1D>();
      }

      // shift and atom index
      const double current_shift( util::ConvertStringToNumericalValue< double>( signal_split( 0)));
      const size_t current_atom_index( util::ConvertStringToNumericalValue< size_t>( signal_split( 2)));

      // acquire the pattern
      // find the first character that is alpha
      const std::string::const_iterator first_pattern_itr( std::find_if( signal_split( 1).begin(), signal_split( 1).end(), std::ptr_fun( &isalpha)));
      const std::string coupling_const_string( signal_split( 1).begin(), first_pattern_itr);
      const std::string pattern_string( first_pattern_itr, signal_split( 1).end());

      const double coupling( util::IsNumerical( coupling_const_string) ? util::ConvertStringToNumericalValue< double>( coupling_const_string) : 0.0);
      const PatternType pattern( pattern_string.empty() ? e_UnknownPattern : PatternTypeFromSingleLetterCode( std::toupper( pattern_string[ pattern_string.length() - 1])));

      static const double s_default_integral( 1.0);

      ATOMS.GotoPosition( current_atom_index);
      // create signal 1D
      return util::ShPtr< Signal1D>
        (
          new Signal1D
          (
            current_shift,
            s_default_integral,
            storage::Vector< storage::Pair< PatternTypeEnum, double> >
            (
              1,
              storage::Pair< PatternTypeEnum, double>( pattern, coupling)
            ),
            *ATOMS
          )
        );
    }

  } // namespace nmr
} // namespace bcl

