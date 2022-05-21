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
#include "score/bcl_score_restraint_xlink.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "restraint/bcl_restraint_handler_atom_distance_assigned.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_restraint_xlink.cpp
  //! @brief this example tests the implementation of the class scoring the agreement of a protein model with data
  //! obtained from a cross-linking experiment
  //!
  //! @author fischea
  //! @date Jun 12, 2014
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreRestraintXlink :
     public ExampleInterface
  {

  //////////
  // data //
  //////////

  public:

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief returns a pointer to a new ExampleScoreRestraintXlink
    //! @return pointer to a new ExampleScoreRestraintXlink
    ExampleScoreRestraintXlink *Clone() const
    {
      return new ExampleScoreRestraintXlink( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief performs the execution of the example
    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // read in cross-linking restraints from an example file
      const std::string restraint_file( AddExampleInputPathToFilename( e_Restraint, "1JL1.xlink_bcl"));
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, restraint_file);
      util::ShPtr< util::ShPtrVector< restraint::AtomDistance> > sp_restraints
      (
        new util::ShPtrVector< restraint::AtomDistance>
        (
          restraint::HandlerAtomDistanceAssigned().ReadRestraints( read)
        )
      );
      io::File::CloseClearFStream( read);

      // construct a scoring object
      const std::string &scheme( "xlink_test");
      const score::RestraintXlink xlink_score( sp_restraints, 5.0, scheme);

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck( xlink_score.GetClassIdentifier(), GetStaticClassName< score::RestraintXlink>());

      // check the scheme of the score
      BCL_ExampleCheck( scheme.compare( xlink_score.GetScheme()), 0);

    ////////////////
    // operations //
    ////////////////

      // create a protein model for testing
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Restraint, "1JL1A.pdb"));
      BCL_ExampleMustOpenInputFile( read, pdb_filename);
      const pdb::Handler pdb( read);
      io::File::CloseClearFStream( read);
      const assemble::ProteinModel model( pdb::Factory().ProteinModelFromPDB( pdb));

      // score the agreement of the protein model with the cross-linking restraints
      const double score_correct_restraints( xlink_score( model));

      // read in incorrect cross-linking restraints from an example file
      const std::string restraint_file_inc( AddExampleInputPathToFilename( e_Restraint, "1JL1_wrong.xlink_bcl"));
      BCL_ExampleMustOpenInputFile( read, restraint_file_inc);
      util::ShPtr< util::ShPtrVector< restraint::AtomDistance> > sp_restraints_incorrect
      (
        new util::ShPtrVector< restraint::AtomDistance>
        (
          restraint::HandlerAtomDistanceAssigned().ReadRestraints( read)
        )
      );
      io::File::CloseClearFStream( read);

      // construct a scoring object for the incorrect restraints
      const score::RestraintXlink xlink_score_incorrect( sp_restraints_incorrect);

      // score the agreement of the protein model with the incorrect cross-linking restraints
      const double score_incorrect_restraints( xlink_score_incorrect( model));

      // make sure that the incorrect restraints score worse than the correct restraints
      BCL_ExampleCheck( score_incorrect_restraints > score_correct_restraints, true);

    //////////////////////
    // input and output //
    //////////////////////

//      // write scoring object
//      WriteBCLObject( xlink_score);
//
//      // read scoring object
//      score::RestraintXlink xlink_score_read;
//      ReadBCLObject( xlink_score_read);
//
//      // compare state and behavior of read in object to written out object
//      BCL_ExampleCheck( scheme.compare( xlink_score_read.GetScheme()), 0);
//      const double score_read( xlink_score_read( model));
//      BCL_ExampleCheck( score_read == score_correct_restraints, true);

      return 0;
    }

  }; // class ExampleScoreRestraintXlink

  //! single instance of this class
  const ExampleClass::EnumType ExampleScoreRestraintXlink::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreRestraintXlink())
  );

} // namespace bcl
