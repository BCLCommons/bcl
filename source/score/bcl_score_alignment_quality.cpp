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
#include "score/bcl_score_alignment_quality.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_quality.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AlignmentQuality::s_Instance
    ( GetObjectInstances().AddInstance( new AlignmentQuality()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AlignmentQuality::AlignmentQuality()
    {
    }

    //! @brief Clone function
    //! @return pointer to new AlignmentQuality
    AlignmentQuality *AlignmentQuality::Clone() const
    {
      return new AlignmentQuality( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AlignmentQuality::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the measure
    //! @param SP_MEASURE the measure to use
    void AlignmentQuality::SetMeasure( const util::ShPtr< quality::MeasureInterface> &SP_MEASURE)
    {
      m_QualityMeasure = SP_MEASURE;
    }

    //! @brief set the atoms to consider for the quality measure
    //! @param ATOM_TYPES storage set of atoms to use
    void AlignmentQuality::SetAtoms( const storage::Set< biol::AtomType> &ATOM_TYPES)
    {
      m_AtomTypes = ATOM_TYPES;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that calculates the quality measure for the given alignment
    //! @param ALIGNMENT the alignment to consider
    //! @return the quality measure for all aligned atom coordinates
    double AlignmentQuality::operator()( const align::AlignmentInterface< biol::AABase> &ALIGNMENT) const
    {
      const storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > coords
      (
        assemble::Quality::CoordinatesFromAlignment( ALIGNMENT, m_AtomTypes)
      );
      return m_QualityMeasure->CalculateMeasure( coords.First(), coords.Second());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AlignmentQuality::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_QualityMeasure, ISTREAM);
      io::Serialize::Read( m_AtomTypes     , ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AlignmentQuality::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_QualityMeasure, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomTypes     , OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
