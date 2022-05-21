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
#include "chemistry/bcl_chemistry_reaction_complete.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"
#include "sdf/bcl_sdf_rxn_factory.h"
#include "sdf/bcl_sdf_rxn_handler.h"

// external includes - sorted alphabetically
#include <cstdio>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_reaction_complete.cpp
  //! @details Tests ReactionComplete class which contains reaction data
  //!
  //! @author geanesar
  //! @date 06/13/2015
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryReactionComplete :
    public ExampleInterface
  {
  public:

    ExampleChemistryReactionComplete *Clone() const
    {
      return new ExampleChemistryReactionComplete( *this);
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

      // read a reaction from file
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "multi_reaction_1.rxn"));
      sdf::RXNHandler handler( input);
      io::File::CloseClearFStream( input);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      chemistry::ReactionComplete new_reaction
      (
        sdf::RXNFactory::MakeReactionComplete( handler)
      );

      util::ShPtr< chemistry::ReactionComplete> sp_rxn
      (
        util::CloneToShPtr( new_reaction)
      );

      BCL_ExampleCheck( new_reaction.GetNumberReactants(), 2);
      BCL_ExampleCheck( new_reaction.GetNumberProducts(), 2);
      BCL_ExampleCheck( sp_rxn->GetNumberReactants(), 2);
      BCL_ExampleCheck( sp_rxn->GetNumberProducts(), 2);

      storage::Map< size_t, storage::Pair< size_t, size_t> > atom_map_r;
      atom_map_r[1] = storage::Pair< size_t, size_t>( 0, 1);
      atom_map_r[2] = storage::Pair< size_t, size_t>( 0, 0);
      atom_map_r[3] = storage::Pair< size_t, size_t>( 0, 2);
      atom_map_r[4] = storage::Pair< size_t, size_t>( 1, 1);

      BCL_ExampleIndirectCheck( new_reaction.GetReactiveAtomsAllReactants(), atom_map_r, "Checking atom mapping was converted correctly for reactants");

      storage::Map< size_t, storage::Pair< size_t, size_t> > atom_map_p;
      atom_map_p[1] = storage::Pair< size_t, size_t>( 0, 0);
      atom_map_p[2] = storage::Pair< size_t, size_t>( 0, 1);
      atom_map_p[3] = storage::Pair< size_t, size_t>( 1, 1);
      atom_map_p[4] = storage::Pair< size_t, size_t>( 0, 2);

      BCL_ExampleIndirectCheck( new_reaction.GetReactiveAtomsAllProducts(), atom_map_p, "Checking atom mapping was converted correctly for products");

      // Equality operator
      BCL_ExampleCheck( new_reaction == *sp_rxn, true);

      // Destroy ReactionComplete
      sp_rxn.Reset();

      // Data access
      std::string output_filename
      ( 
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "reaction_complete.rxn")
      );
      std::string symmetry_filename( output_filename + ".symmetry");
      std::string correct_filename( output_filename + ".correct");
      
      // Test the correctness of writing the reaction
      io::OFStream output;
      BCL_ExampleMustOpenOutputFile( output, output_filename);

      // test read/write symmetry
      new_reaction.WriteRXN( output);
      io::File::CloseClearFStream( output);

      // Check symmetry of read/write by reading in the just-written file
      BCL_ExampleMustOpenInputFile( input, output_filename);
      sdf::RXNHandler sym_handler( input);
      chemistry::ReactionComplete second_rxn( sdf::RXNFactory::MakeReactionComplete( sym_handler));
      io::File::CloseClearFStream( input);

      BCL_ExampleMustOpenOutputFile( output, symmetry_filename)
      second_rxn.WriteRXN( output);
      io::File::CloseClearFStream( output);

      // check the symmetry file matches the original output file
      if( BCL_ExampleCheck( io::File::FilesMatch( output_filename, symmetry_filename), true))
      {
        remove( symmetry_filename.c_str());
      }
      else
      {
        BCL_MessageStd
        ( 
          "File output by example, \"" + output_filename + "\", did not match target file \"" + symmetry_filename + "\""
        );
      }
      
      // check if the output file is right
      if( BCL_ExampleCheck( io::File::FilesMatch( output_filename, correct_filename), true)
      )
      {
        remove( output_filename.c_str());
      }
      else
      {
        BCL_MessageStd
        ( 
          "File output by example, \"" + output_filename + "\", did not match target file \"" + correct_filename + "\""
        );
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSmallMolecule

  const ExampleClass::EnumType ExampleChemistryReactionComplete::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryReactionComplete())
  );

} // namespace bcl
