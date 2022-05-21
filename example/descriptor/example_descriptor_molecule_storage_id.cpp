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
#include "descriptor/bcl_descriptor_molecule_storage_id.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_molecule_ensemble.h"
#include "chemistry/bcl_chemistry_molecule_storage_file.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_storage_id.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 18, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeStorageId :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeStorageId *Clone() const
    {
      return new ExampleDescriptorMoleculeStorageId( *this);
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
      // create input stream for reading a smallmolecule ensemble
      io::IFStream input;
      const std::string test_file( AddExampleInputPathToFilename( e_Chemistry, "test_set.5_structures.sdf"));
      BCL_ExampleMustOpenInputFile( input, test_file);
      // read in ensemble
      chemistry::MoleculeEnsemble ensemble( input);
      // close stream
      io::File::CloseClearFStream( input);

      descriptor::MoleculeStorageId id_defualt( new chemistry::MoleculeStorageFile( true));

      // make a new storage id calculator via a small molecule property
      const std::string label( id_defualt.GetAlias() + "(" + test_file + ")");
      descriptor::CheminfoID id_calc( label);

      // check that the molecule property existed
      BCL_ExampleIndirectCheck
      (
        descriptor::CheminfoID( label).IsDefined(),
        true,
        "MoleculeStorageId construction from: " + label
      );
      BCL_ExampleIndirectCheck
      (
        descriptor::CheminfoID( label)->GetSizeOfFeatures(),
        14,
        "MoleculeStorageId construction from: " + label
      );

      // check that the molecule ids are calculated properly
      size_t index( 0);
      bool keys_were_right( true);
      for
      (
        storage::List< chemistry::MoleculeComplete>::const_iterator itr( ensemble.Begin()), itr_end( ensemble.End());
        itr != itr_end && keys_were_right;
        ++itr, ++index
      )
      {
        std::string id_vector( util::TrimString( std::string( id_calc->SumOverObject( *itr).Begin(), 14)));
        if( id_vector != util::Format()( index))
        {
          // only throw an example check if the molecule id was wrong; to avoid artificially inflating the number
          // of example checks reported for this example
          BCL_ExampleIndirectCheck
          (
            id_vector,
            util::Format()( index),
            "MoleculeStorageId determined keys properly from a file"
          );
          keys_were_right = false;
        }
      }
      // only perform the example check that the keys were right if they were right; otherwise the BCL_ExampleCheck
      // in the for loop caught it
      if( keys_were_right)
      {
        BCL_ExampleIndirectCheck
        (
          keys_were_right,
          true,
          "MoleculeStorageId determined keys properly from a file"
        );
      }

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeStorageId

  const ExampleClass::EnumType ExampleDescriptorMoleculeStorageId::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeStorageId())
  );

} // namespace bcl
