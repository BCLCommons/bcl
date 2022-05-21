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
#include "score/bcl_score_restraint_atom_distance.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "nmr/bcl_nmr_star_noe_handler.h"
#include "restraint/bcl_restraint_atom_distance.h"
#include "score/bcl_score_restraint_noe_knowledge_based.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_restraint_atom_distance.cpp
  //!
  //! @author weinerbe
  //! @date Mar 16, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreRestraintAtomDistance :
    public ExampleInterface
  {
  public:

    ExampleScoreRestraintAtomDistance *Clone() const
    {
      return new ExampleScoreRestraintAtomDistance( *this);
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
      // create strings for input files
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "2KDC_single_model.pdb"));
      const std::string restraint_file( AddExampleInputPathToFilename( e_Biology, "nmr_star_restraints.txt"));

      // create protein model from pdb
      BCL_MessageStd( "building model");
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone)
      );

      // read restraints
      nmr::StarNOEHandler handler;
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, restraint_file);

      // create restraints
      util::ShPtr< assemble::ProteinModelData> sp_model_data( protein_model.GetProteinModelData());
      util::ShPtrVector< restraint::AtomDistance> atom_distances( handler.ReadRestraints( read));

      io::File::CloseClearFStream( read);
      protein_model.SetProteinModelData( sp_model_data);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test GetDefaultScheme
      const std::string def_scheme( "distance_restraint");
      BCL_ExampleCheck( score::RestraintAtomDistance::GetDefaultScheme(), def_scheme);

      // test default constructor
      score::RestraintAtomDistance def_construct;
      BCL_ExampleIndirectCheck( def_construct.GetScheme(), def_scheme, "default constructor");

      // test constructor from a scoring function and a scheme
      score::RestraintAtomDistance score_construct
      (
        score::RestraintNoeKnowledgeBased(),
        1.0,
        "test_noe_score"
      );
      score_construct.SetRestraints( util::CloneToShPtr( atom_distances));

    /////////////////
    // data access //
    /////////////////

      // test GetScheme
      BCL_ExampleCheck( score_construct.GetScheme(), "test_noe_score");

    ///////////////
    // operators //
    ///////////////

      // test () operator
      const double correct_score( -43.8393);
      const double calculated_score( score_construct( protein_model));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( calculated_score, correct_score),
        true,
        "() operator should return " + util::Format()( correct_score) +
        " but instead returns " + util::Format()( calculated_score)
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write out the detailed scheme and values
      score_construct.WriteDetailedSchemeAndValues( protein_model, util::GetLogger());

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreRestraintAtomDistance

  const ExampleClass::EnumType ExampleScoreRestraintAtomDistance::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreRestraintAtomDistance())
  );

} // namespace bcl
