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
#include "fold/bcl_fold_mutate_sheet_fit_to_template.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_fold_template.h"
#include "math/bcl_math_mutate_result.h"
#include "quality/bcl_quality_rmsd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateSheetFitToTemplate::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateSheetFitToTemplate())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief returns a pointer to a new MutateSheetFitToTemplate
    //! @return pointer to a new MutateSheetFitToTemplate
    MutateSheetFitToTemplate *MutateSheetFitToTemplate::Clone() const
    {
      return new MutateSheetFitToTemplate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &MutateSheetFitToTemplate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief fits the given sheet to a randomly selected template
    //! @param SHEET sheet to be fitted
    //! @return the fitted sheet and a boolean indicating if the application of the mutate was successful
    math::MutateResult< assemble::Domain> MutateSheetFitToTemplate::operator()( const assemble::Domain &SHEET) const
    {
      BCL_MessageTop( "Fitting sheet to template");

      // initialize a undefined domain
      static util::ShPtr< assemble::Domain> s_sp_undefined_domain;

      // make sure the passed domain is a sheet and has at least 2 strands
      if
      (
        SHEET.GetNumberSSEs() < 2 &&
        (
          SHEET.GetTopology()->GetType() != assemble::Topology::e_Sheet ||
          SHEET.GetTopology()->GetType() != assemble::Topology::e_BetaBarrel
        )
      )
      {
        BCL_MessageVrb( "The provided sheet has less than 2 strands or is not a sheet, skipping.");
        return math::MutateResult< assemble::Domain>( s_sp_undefined_domain, *this);
      }

      // get the ordered vector of strands in the sheet and cast to SSE pointers from SSEGeometryInterface pointers
      util::SiPtrVector< const assemble::SSE> ordered_strands( SHEET.GetTopology()->GetElements());

      // make sure the cast was successful
      if( !ordered_strands.IsDefined())
      {
        BCL_MessageCrt( "The cast to SSE pointers from SSEGeometryInterface pointers failed!");
        return math::MutateResult< assemble::Domain>( s_sp_undefined_domain, *this);
      }

      // randomly reverse the vector
      if( random::GetGlobalRandom().Boolean())
      {
        util::SiPtrVector< const assemble::SSE> temp_vector
        (
          ordered_strands.ReverseBegin(),
          ordered_strands.ReverseEnd()
        );
        ordered_strands = temp_vector;
      }

      // get a random template that has the same number of strands
      const assemble::FoldTemplate &fold_template( assemble::SheetTemplateHandler::GetRandomTemplate( ordered_strands));

      // application of the mutate failed if no suitable template was found
      if( fold_template.GetGeometries().IsEmpty())
      {
        return math::MutateResult< assemble::Domain>( s_sp_undefined_domain, *this);
      }

      // fit the SSEs to the template
      util::ShPtr< assemble::Domain> sp_new_sheet( fold_template.FitSSEs( ordered_strands).Clone());

      // now we need to move the fitted strands close to their original locations
      // first calculate the transformation needed
      math::TransformationMatrix3D transformation
      (
        quality::RMSD::SuperimposeCoordinates
        (
          SHEET.GetAtomCoordinates(),
          sp_new_sheet->GetAtomCoordinates()
        )
      );

      // make sure the transformation is defined
      if( !transformation.IsDefined())
      {
        BCL_MessageCrt( "The transformation of the sheet failed!");
        return math::MutateResult< assemble::Domain>( s_sp_undefined_domain, *this);
      }

      // transform the new strands
      sp_new_sheet->Transform( transformation);

      return math::MutateResult< assemble::Domain>( sp_new_sheet, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read object from input stream
    //! @param ISTREAM input stream to read object from
    //! @return input stream which was read from
    std::istream &MutateSheetFitToTemplate::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write object into  output stream
    //! @param OSTREAM output stream to write object into
    //! @param INDENT number of indentations to separate members
    //! @return output stream object was written into
    std::ostream &MutateSheetFitToTemplate::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
