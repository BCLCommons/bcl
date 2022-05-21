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

#ifndef BCL_SCORE_ALIGNMENT_QUALITY_H_
#define BCL_SCORE_ALIGNMENT_QUALITY_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "align/bcl_align.fwd.hh"
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom_types.h"
#include "function/bcl_function_unary_interface.h"
#include "quality/bcl_quality_measure_interface.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignmentQuality
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_score_alignment_quality.cpp @endlink
    //! @author woetzen
    //! @date Mar 14, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AlignmentQuality :
      public function::UnaryInterface< const align::AlignmentInterface< biol::AABase>, double>
    {

    private:

    //////////
    // data //
    //////////

      //! quality measure to use
      util::ShPtr< quality::MeasureInterface> m_QualityMeasure;

      //! atoms to use
      storage::Set< biol::AtomType> m_AtomTypes;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AlignmentQuality();

      //! @brief Clone function
      //! @return pointer to new AlignmentQuality
      AlignmentQuality *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief set the measure
      //! @param SP_MEASURE the measure to use
      void SetMeasure( const util::ShPtr< quality::MeasureInterface> &SP_MEASURE);

      //! @brief set the atoms to consider for the quality measure
      //! @param ATOM_TYPES storage set of atoms to use
      void SetAtoms( const storage::Set< biol::AtomType> &ATOM_TYPES);

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that calculates the quality measure for the given alignment
      //! @param ALIGNMENT the alignment to consider
      //! @return the quality measure for all aligned atom coordinates
      double operator()( const align::AlignmentInterface< biol::AABase> &ALIGNMENT) const;

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

    }; // class AlignmentQuality

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_ALIGNMENT_QUALITY_H_
