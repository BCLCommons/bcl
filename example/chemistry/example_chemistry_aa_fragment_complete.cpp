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
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_aa_fragment_complete.cpp
  //!
  //! @author mendenjl
  //! @date Sep 16, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryAAFragmentComplete :
    public ExampleInterface
  {
  public:

    ExampleChemistryAAFragmentComplete *Clone() const
    {
      return new ExampleChemistryAAFragmentComplete( *this);
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

      // initialize sdf filename
      pdb::Factory factory( biol::GetAAClasses().e_AAComplete);
      assemble::ProteinModel model
      (
        factory.ProteinModelFromPDBFilename
        (
          AddExampleInputPathToFilename( e_Biology, "2HAC.pdb"),
          storage::Map< biol::SSType, size_t>::Create
          (
            std::make_pair( biol::GetSSTypes().HELIX, size_t( 0)),
            std::make_pair( biol::GetSSTypes().STRAND, size_t( 0)),
            std::make_pair( biol::GetSSTypes().COIL, size_t( 0))
          )
        )
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      chemistry::AAFragmentComplete aa_fragment( model.GetAminoAcids());

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      aa_fragment.StoreProperty
      (
        descriptor::GetCheminfoProperties().calc_SigmaCharge.GetString(),
        descriptor::GetCheminfoProperties().calc_SigmaCharge->CollectValuesOnEachElementOfObject( aa_fragment)
      );
      aa_fragment.StoreProperty
      (
        descriptor::GetCheminfoProperties().calc_PiCharge.GetString(),
        descriptor::GetCheminfoProperties().calc_PiCharge->CollectValuesOnEachElementOfObject( aa_fragment)
      );
      aa_fragment.StoreProperty
      (
        descriptor::GetCheminfoProperties().calc_VCharge.GetString(),
        descriptor::GetCheminfoProperties().calc_VCharge->CollectValuesOnEachElementOfObject( aa_fragment)
      );

      BCL_ExampleCheck( aa_fragment.GetSize(), 550);

      // uncomment this line to print out the molecule
      //aa_fragment.WriteMDL( util::GetLogger());
      pdb::Factory().WriteModelToPDB( aa_fragment.ReconstructProteinModel( model), util::GetLogger());
      assemble::ProteinModel model2
      (
        factory.ProteinModelFromPDBFilename
        (
          AddExampleInputPathToFilename( e_Biology, "2HAC.pdb"),
          storage::Map< biol::SSType, size_t>::Create
          (
            std::make_pair( biol::GetSSTypes().HELIX, size_t( 0)),
            std::make_pair( biol::GetSSTypes().STRAND, size_t( 0)),
            std::make_pair( biol::GetSSTypes().COIL, size_t( 999))
          )
        )
      );
      pdb::Factory().WriteModelToPDB( chemistry::AAFragmentComplete( model2.GetAminoAcids(), true).ReconstructProteinModel( model2), util::GetLogger());

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryAAFragmentComplete

  const ExampleClass::EnumType ExampleChemistryAAFragmentComplete::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryAAFragmentComplete())
  );

} // namespace bcl
