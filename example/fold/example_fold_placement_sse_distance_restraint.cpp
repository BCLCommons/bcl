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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "fold/bcl_fold_placement_sse_distance_restraint.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "nmr/bcl_nmr_star_noe_handler.h"
#include "restraint/bcl_restraint_atom_distance.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_placement_sse_distance_restraint.cpp
  //!
  //! @author weinerbe
  //! @date Jan 11, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldPlacementSSEDistanceRestraint :
    public ExampleInterface
  {
  public:

    ExampleFoldPlacementSSEDistanceRestraint *Clone() const
    {
      return new ExampleFoldPlacementSSEDistanceRestraint( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // read in the protein model
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel
        (
          AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"),
          biol::GetAAClasses().e_AABackBone,
          ssetype_min_size
        )
      );

      // read in the restraints
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1ubi.noe_star"));
      nmr::StarNOEHandler handler;
      util::ShPtrVector< restraint::AtomDistance> atom_distances( handler.ReadRestraints( read));
      io::File::CloseClearFStream( read);

      // get the helix
      util::ShPtr< assemble::SSE> sp_helix
      (
        protein_model.GetSSEs( biol::GetSSTypes().HELIX).FirstElement()->Clone()
      );
      sp_helix->SetToIdealConformationAtOrigin();
      protein_model.Remove( *sp_helix);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      fold::PlacementSSEDistanceRestraint def_construct;

      // test constructor from NOE type
      const fold::PlacementSSEDistanceRestraint noe_construct( util::CloneToShPtr( atom_distances));

    ////////////////
    // operations //
    ////////////////

      // test Place
      BCL_ExampleCheck( def_construct.Place( *sp_helix, protein_model).Second(), false);
      BCL_ExampleCheck( noe_construct.Place( *sp_helix, protein_model).Second(), true);

    //////////////////////
    // input and output //
    //////////////////////

      // test write and read
      WriteBCLObject( noe_construct);
      ReadBCLObject( def_construct);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldPlacementSSEDistanceRestraint

  const ExampleClass::EnumType ExampleFoldPlacementSSEDistanceRestraint::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldPlacementSSEDistanceRestraint())
  );

} // namespace bcl
