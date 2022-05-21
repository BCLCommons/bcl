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
#include "chemistry/bcl_chemistry_reaction_structure.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_types.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_atom_info.h"
#include "sdf/bcl_sdf_fragment_factory.h"
#include "sdf/bcl_sdf_rxn_factory.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_reaction_structure.cpp
  //! @details Tests ChemistryReactionStructure class which contains small molecule configuration data
  //!
  //! @author geanesar
  //! @date
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryReactionStructure :
    public ExampleInterface
  {
  public:

    ExampleChemistryReactionStructure *Clone() const
    {
      return new ExampleChemistryReactionStructure( *this);
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
      std::string template_filename
      ( 
        AddExampleInputPathToFilename
        ( 
          e_Chemistry, "reaction_structure_template.sdf"
        )
      );

      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, template_filename);

      sdf::MdlHandler handler1( input);
      BCL_ExampleIndirectAssert( handler1.IsValid(), true, "Requisite MdlHandler 1 validity");
      sdf::MdlHandler handler2( input);
      BCL_ExampleIndirectAssert( handler2.IsValid(), true, "Requisite MdlHandler 2 validity");

      io::File::CloseClearFStream( input);

      chemistry::ReactionStructure arom_template( sdf::RXNFactory::MakeReactionStructure( handler1.GetMolfile()));
      chemistry::ReactionStructure aliph_template( sdf::RXNFactory::MakeReactionStructure( handler2.GetMolfile()));

      std::string query_filename
      ( 
        AddExampleInputPathToFilename
        ( 
          e_Chemistry, "reaction_structure_test.sdf"
        )
      );
      
      BCL_ExampleMustOpenInputFile( input, query_filename);

      sdf::MdlHandler test_handler1( input);
      BCL_ExampleIndirectAssert( test_handler1.IsValid(), true, "Requisite MdlHandler 3 validity");
      sdf::MdlHandler test_handler2( input);
      BCL_ExampleIndirectAssert( test_handler2.IsValid(), true, "Requisite MdlHandler 4 validity");

      chemistry::FragmentComplete arom_mol( sdf::FragmentFactory::MakeFragment( test_handler1));
      chemistry::FragmentComplete aliph_mol( sdf::FragmentFactory::MakeFragment( test_handler2));

      io::File::CloseClearFStream( input);

      BCL_ExampleCheck( arom_template.GetSize(), 2);
      BCL_ExampleCheck( aliph_template.GetSize(), 2);
      
      // aromatic molecule (mol 1) matches the aromatic template (template 1)
      BCL_ExampleCheck( arom_template.ContainedIn( arom_mol), true);

      // aliphatic molecule (mol 2) matches the aliphatic template (template 2)
      BCL_ExampleCheck( aliph_template.ContainedIn( aliph_mol), true);

      // aromatic molecule (mol 1) does not match the aliphatic template (template 2)
      BCL_ExampleCheck( aliph_template.ContainedIn( arom_mol), false);

      // aromatic molecule (mol 2) does not match the aliphatic template (template 1)
      BCL_ExampleCheck( arom_template.ContainedIn( aliph_mol), false);
      
      return 0;
    }
    // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSmallMolecule

  const ExampleClass::EnumType ExampleChemistryReactionStructure::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryReactionStructure())
  );

} // namespace bcl
