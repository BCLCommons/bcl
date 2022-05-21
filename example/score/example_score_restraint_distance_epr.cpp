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
#include "score/bcl_score_restraint_distance_epr.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "restraint/bcl_restraint_handler_atom_distance_assigned.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_restraint_distance_epr.cpp
  //! @brief this example tests the implementation of the class scoring the agreement of a protein model with distance
  //! restraints obtained from an EPR experiment
  //!
  //! @author fischea
  //! @date Jun 7, 2014
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreRestraintDistanceEPR :
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

    //! @brief returns a pointer to a new ExampleScoreRestraintDistanceEPR
    //! @return pointer to a new ExampleScoreRestraintDistanceEPR
    ExampleScoreRestraintDistanceEPR *Clone() const
    {
      return new ExampleScoreRestraintDistanceEPR( *this);
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

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // read in EPR restraints from an example file
      const std::string restraint_file( AddExampleInputPathToFilename( e_Restraint, "1F16.epr_cst_bcl"));
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, restraint_file);
      util::ShPtrVector< restraint::AtomDistance> restraints
      (
        restraint::HandlerAtomDistanceAssigned().ReadRestraints( read)
      );
      io::File::CloseClearFStream( read);

      // construct scoring object from restraints, histogram, and scheme
      const std::string histogram_file( "sl-cb_distances_geom.histogram");
      const std::string scheme( "epr_scheme");
      score::RestraintDistanceEPR epr_score( restraints, histogram_file, scheme);

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck( epr_score.GetClassIdentifier(), GetStaticClassName< score::RestraintDistanceEPR>());

      // check the scheme of the score
      BCL_ExampleCheck( scheme.compare( epr_score.GetScheme()), 0);

    ////////////////
    // operations //
    ////////////////

      // create a protein model for testing
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Restraint, "1F16_8.pdb"));
      BCL_ExampleMustOpenInputFile( read, pdb_filename);
      pdb::Handler pdb( read);
      io::File::CloseClearFStream( read);
      assemble::ProteinModel model( pdb::Factory().ProteinModelFromPDB( pdb));

      // score the agreement of the protein model with the EPR restraints
//      const double score( epr_score( model));

      return 0;
    }

  }; // class ExampleScoreRestraintDistanceEPR

  //! single instance of this class
  const ExampleClass::EnumType ExampleScoreRestraintDistanceEPR::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreRestraintDistanceEPR())
  );

} // namespace bcl
