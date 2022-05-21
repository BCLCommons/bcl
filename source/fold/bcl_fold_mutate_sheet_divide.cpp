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
#include "fold/bcl_fold_mutate_sheet_divide.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateSheetDivide::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateSheetDivide( 4, 2, false))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from min sheet and divided sheet sizes, and minimum and maximum translations
    //! @param MIN_SHEET_SIZE minimum number of strands in the sheet to be divided
    //! @param MIN_DIVIDED_SHEET_SIZE minimum number of strand in the divided sheet
    //! @param FORM_BETA_SANDWICH whether to form a beta sandwich
    //! @param MIN_TRANSLATIONS minimum translation along x,y and z axis
    //! @param MAX_TRANSLATIONS maximum translation along x,y and z axis
    MutateSheetDivide::MutateSheetDivide
    (
      const size_t MIN_SHEET_SIZE,
      const size_t MIN_DIVIDED_SHEET_SIZE,
      const bool FORM_BETA_SANDWICH,
      const linal::Vector3D &MIN_TRANSLATIONS,
      const linal::Vector3D &MAX_TRANSLATIONS
    ) :
      m_MinSheetSize( MIN_SHEET_SIZE),
      m_MinDividedSheetSize( MIN_DIVIDED_SHEET_SIZE),
      m_FormBetaSandwich( FORM_BETA_SANDWICH),
      m_MinTranslations( MIN_TRANSLATIONS),
      m_MaxTranslations( MAX_TRANSLATIONS)
    {
      // make sure the given sheet size large enough to make sure both sheets after division are at least m_MinDividedSheetSize
      BCL_Assert
      (
        MIN_SHEET_SIZE >= 2 * MIN_DIVIDED_SHEET_SIZE,
        "The given min sheet size " + util::Format()( MIN_SHEET_SIZE) +
        " should be at least twice the min divided sheet size " + util::Format()( MIN_DIVIDED_SHEET_SIZE)
      );

      // make sure the min divided sheet size is at least two
      BCL_Assert
      (
        MIN_DIVIDED_SHEET_SIZE >= 1,
        "The given min divided sheet size should be at least 1 not " + util::Format()( MIN_DIVIDED_SHEET_SIZE)
      );
    }

    //! @brief Clone function
    //! @return pointer to new MutateSheetDivide
    MutateSheetDivide *MutateSheetDivide::Clone() const
    {
      return new MutateSheetDivide( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateSheetDivide::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that takes a Sheet and return a mutated Sheet
    //! @param SHEET Sheet which will be mutated
    //! @return MutateResult with the mutated Sheet
    math::MutateResult< assemble::Domain> MutateSheetDivide::operator()( const assemble::Domain &SHEET) const
    {
      // initialize empty sheet to return
      static util::ShPtr< assemble::Domain> s_empty_sheet;

      // make sure the passed domain has a valid topology and is of type sheet or beta-barrel
      if
      (
        !SHEET.GetTopology().IsDefined() ||
        !(
           SHEET.GetTopology()->GetType() == assemble::Topology::e_Sheet ||
           SHEET.GetTopology()->GetType() == assemble::Topology::e_BetaBarrel
         )
      )
      {
        // warn user and return
        BCL_MessageVrb( "The given domain is not a sheet or a barrel");
        return math::MutateResult< assemble::Domain>( s_empty_sheet, *this);
      }

      // store the total number of strand
      const size_t total_nr_strands( SHEET.GetNumberSSEs());

      // if not enough strands
      if( total_nr_strands < m_MinSheetSize)
      {
        BCL_MessageDbg( "The sheet has less than minimum required number of strands");
        return math::MutateResult< assemble::Domain>( s_empty_sheet, *this);
      }

      // create a new sheet
      util::ShPtr< assemble::Domain> new_sheet( SHEET.Clone());

      // make a copy of the order vector
      util::SiPtrVector< const assemble::SSEGeometryInterface> order_vector( SHEET.GetTopology()->GetElements());

      // locate the dividing index
      const size_t dividing_index
      (
        random::GetGlobalRandom().SizeT
        (
          math::Range< size_t>( m_MinDividedSheetSize, total_nr_strands - m_MinDividedSheetSize)
        )
      );

      // random boolean left or right
      const bool left_or_right( random::GetGlobalRandom().Boolean());

      // construct the new vector
      util::SiPtrVector< const assemble::SSEGeometryInterface> subsheet
      (
        left_or_right ?
          order_vector.SubSiPtrVector( 0, dividing_index) :
          order_vector.SubSiPtrVector( dividing_index, total_nr_strands - dividing_index)
      );

      // construct the vector that contains the part of the sheet that will stay intact
      util::SiPtrVector< const assemble::SSEGeometryInterface> subsheet_intact
      (
        left_or_right ?
          order_vector.SubSiPtrVector( dividing_index, total_nr_strands - dividing_index) :
          order_vector.SubSiPtrVector( 0, dividing_index)
      );

      // initialize translation that will be applied
      linal::Vector3D translation;

      // if sandwich
      if( m_FormBetaSandwich)
      {
        // calculate the center of the subsheet intact and the subsheet
        linal::Vector3D center_subsheet_intact( subsheet_intact( ( subsheet_intact.GetSize() - 1) / 2)->GetCenter());
        linal::Vector3D center_subsheet( subsheet( ( subsheet.GetSize() - 1) / 2)->GetCenter());

        // calculate the translation that will move the center of subsheet to center of subsheet_intact
        translation += ( center_subsheet_intact - center_subsheet);

        // now determine the x translation that will be move sheet apart
        // this can be determined randomly in the of preferred distances
        const linal::Vector3D x_translation
        (
          SHEET.GetOrientation().GetAxis( coord::GetAxes().e_X) *
          random::GetGlobalRandom().Sign() *
          random::GetGlobalRandom().Double
          (
            contact::GetTypes().SHEET_SHEET->GetPreferredDistanceRange()
          )
        );

        // sum up the translations
        translation += x_translation;
      }
      // else
      else
      {

        // determine the order of y axis
        // this needs to point away from the sheet so that y translation does not overlap the subsheet with
        // the other part of the sheet that is not selected
        // the y axis of a sheet by default points from left to right in the sheet order vector
        // therefore if we picked the left side
        // then we want to translate in the opposite direction of the y axis of the sheet
        const linal::Vector3D y_axis
        (
          left_or_right ?
            -SHEET.GetOrientation().GetAxis( coord::GetAxes().e_Y) :
            SHEET.GetOrientation().GetAxis( coord::GetAxes().e_Y)
        );

        // now determine the translations in each axis
        // x and z translation can be in either direction
        const linal::Vector3D x_translation
        (
          SHEET.GetOrientation().GetAxis( coord::GetAxes().e_X) *
          random::GetGlobalRandom().Sign() *
          random::GetGlobalRandom().Double
          (
            math::Range< double>( m_MinTranslations( coord::GetAxes().e_X), m_MaxTranslations( coord::GetAxes().e_X))
          )
        );
        const linal::Vector3D z_translation
        (
          SHEET.GetOrientation().GetAxis( coord::GetAxes().e_Z) *
          random::GetGlobalRandom().Sign() *
          random::GetGlobalRandom().Double
          (
            math::Range< double>( m_MinTranslations( coord::GetAxes().e_Z), m_MaxTranslations( coord::GetAxes().e_Z))
          )
        );

        // y translation only occurs in the direction of the calculated axis as explained previously
        const linal::Vector3D y_translation
        (
          SHEET.GetOrientation().GetAxis( coord::GetAxes().e_Y) *
          random::GetGlobalRandom().Double
          (
            math::Range< double>( m_MinTranslations( coord::GetAxes().e_Y), m_MaxTranslations( coord::GetAxes().e_Y))
          )
        );

        // sum up the translations
        translation = ( x_translation + y_translation + z_translation);
      }

      // now iterate over the sses in the subset
      for
      (
        util::SiPtrVector< const assemble::SSEGeometryInterface>::const_iterator
          sse_itr( subsheet.Begin()), sse_itr_end( subsheet.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // make a copy of this sse
        util::ShPtr< assemble::SSE> new_sse( ( *sse_itr)->Clone());
        // check the dynamic cast succeeded
        BCL_Assert( new_sse.IsDefined(), "The dynamic cast from SSEGeometryInterface to SSE failed!");

        // apply the translation
        new_sse->Translate( translation);

        // replace in the new sheet
        new_sheet->Replace( new_sse);
      }

      // end
      return math::MutateResult< assemble::Domain>( new_sheet, *this);

    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateSheetDivide::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_MinSheetSize, ISTREAM);
      io::Serialize::Read( m_MinDividedSheetSize, ISTREAM);
      io::Serialize::Read( m_FormBetaSandwich, ISTREAM);
      io::Serialize::Read( m_MinTranslations, ISTREAM);
      io::Serialize::Read( m_MaxTranslations, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateSheetDivide::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_MinSheetSize, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinDividedSheetSize, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FormBetaSandwich, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinTranslations, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MaxTranslations, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
