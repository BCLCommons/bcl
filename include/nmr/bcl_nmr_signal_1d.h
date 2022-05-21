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

#ifndef BCL_NMR_SIGNAL_1D_H_
#define BCL_NMR_SIGNAL_1D_H_

// include the namespace header
#include "bcl_nmr.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "iterate/bcl_iterate_generic.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Signal1D
    //! @brief This class is used to store the signals of 1 atom, spectra consisting of 1D and ND
    //! signals are constructed using Signal which contains multiple Signal1D.
    //!
    //! @see @link example_nmr_signal_1d.cpp @endlink
    //! @author mueller, butkiem1, woetzen
    //! @date 04/11/2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Signal1D :
      public util::ObjectInterface
    {

    public:

      //! type of 1 dimensional patterns in NMR spectrum
      enum PatternType
      {
        e_UnknownPattern,
        e_Singlet,
        e_Doublet,
        e_Triplet,
        e_Quadruplet,
        e_Multiplet,
        s_NumberPatternTypes
      };

      //! @brief PatternType as string
      //! @param PATTERN_TYPE the PatternType
      //! @return the string for the PatternType
      static const std::string &GetPatternTypeDescriptor( const PatternType &PATTERN_TYPE);

      //! @brief PatternTypeEnum is used for I/O of PatternType
      typedef util::WrapperEnum< PatternType, &GetPatternTypeDescriptor, s_NumberPatternTypes> PatternTypeEnum;

      //! @brief PatternType as single character
      //! @details used in mdl Spectrum lines to describe the pattern
      //! @param PATTERN_TYPE the PatternType
      //! @return the char for the PatternType
      static const char &GetPatternTypeSingleLetterCode( const PatternType &PATTERN_TYPE);

      //! @brief PatternType from single letter code
      //! @details convert char into a PatternType
      //! @param SINGLE_LETTER_CODE the single letter code as it may appear in an mdl Spectrum line
      //! @return PatternType for single letter code, e_UnknownPattern if not found
      static PatternType PatternTypeFromSingleLetterCode( const char SINGLE_LETTER_CODE);

    //////////
    // data //
    //////////

      //! @brief separator for mdl signal line single entry
      static const char s_InformationSeparator = ';';

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      double                                                     m_ChemicalShift; //!< chemical shift in ppm
      double                                                     m_Integral;      //!< integral over the shift
      storage::Vector< storage::Pair< PatternTypeEnum, double> > m_Pattern;       //!< describes the splitting pattern
      util::SiPtr< const chemistry::AtomConformationalInterface> m_Atom;          //!< atom to which the chemical shift belongs

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief standard constructor
      Signal1D()
      {
      }

      //! @brief construct signal1D from its properties
      Signal1D
      (
        const double CHEMICAL_SHIFT,
        const double INTEGRAL,
        const storage::Vector< storage::Pair< PatternTypeEnum, double> > &PATTERN,
        const chemistry::AtomConformationalInterface &ATOM
      ) :
        m_ChemicalShift( CHEMICAL_SHIFT),
        m_Integral( INTEGRAL),
        m_Pattern( PATTERN),
        m_Atom( &ATOM)
      {
      }

      //! @brief virtual copy constructor
      Signal1D *Clone() const
      {
        return new Signal1D( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get chemical shift
      //! @return chemical shift
      double GetChemicalShift() const
      {
        return m_ChemicalShift;
      }

      //! @brief set chemical shift to given CHEMICALSHIFT
      void SetChemicalShift( const double CHEMICAL_SHIFT)
      {
        m_ChemicalShift = CHEMICAL_SHIFT;
      }

      //! @brief get integral
      //! @return Integral of Signal
      double GetIntegral() const
      {
        return m_Integral;
      }

      //! @brief set integral to given INTEGRAL
      void SetIntegral( double const INTEGRAL)
      {
        m_Integral = INTEGRAL;
      }

      //! @brief get atom involved that is involved in signal
      //! @return assigned Atom
      const util::SiPtr< const chemistry::AtomConformationalInterface> &GetAtomInvolvedInSignal() const
      {
        return m_Atom;
      }

      //! @brief set assigned atom to given ATOM
      void SetAtomInvolvedInSignal( const chemistry::AtomConformationalInterface &ATOM)
      {
        m_Atom = &ATOM;
      }

      //! @brief get Signal1D pattern
      //! @return the pattern
      const storage::Vector< storage::Pair< PatternTypeEnum, double> > &GetPattern() const
      {
        return m_Pattern;
      }

      //! @brief set the pattern to given PATTERN
      //! @param PATTERN Signal1D pattern
      void SetPattern( const storage::Vector< storage::Pair< PatternTypeEnum, double> > &PATTERN)
      {
        m_Pattern = PATTERN;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read Spectrum from std::istream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write Spectrum to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief Signal from single entry in mdl
      //! @param SIGNAL_STRING string of the form shift;coupling_constcoupling;atom_index
      //! @param ATOMS the atoms, the atom index is relative to
      //! @return Signal1D for that string
      static util::ShPtr< Signal1D> SignalFromString
      (
        const std::string &SIGNAL_STRING,
        iterate::Generic< const chemistry::AtomConformationalInterface> ATOMS
      );

    }; // class Signal1D

  } // namespace nmr
} // namespace bcl

#endif //BCL_NMR_SIGNAL_1D_H_

