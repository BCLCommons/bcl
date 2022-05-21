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
#include "assemble/bcl_assemble_collector_aa_specified.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_collector_aa_specified.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleCollectorAASpecified :
    public ExampleInterface
  {
  public:

    ExampleAssembleCollectorAASpecified *Clone() const
    {
      return new ExampleAssembleCollectorAASpecified( *this);
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
      // get protein model
      assemble::ProteinModel model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2LZM.pdb")));

      // create string with the filename containing a list of information for creating aa locators
      std::string resi_filename( AddExampleInputPathToFilename( e_Biology, "collector_aa_specified.ls"));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      {
        assemble::CollectorAASpecified def_constr;
        const util::SiPtrList< const biol::AABase> collected_resis( def_constr.Collect( model));
        BCL_ExampleCheck( collected_resis.IsEmpty(), true);
      }

      // constructor taking filename parameter
      assemble::CollectorAASpecified param_constr( resi_filename);
      {
        const util::SiPtrList< const biol::AABase> collected_resis( param_constr.Collect( model));
        BCL_ExampleCheck( collected_resis.GetSize(), 2);
      }

      // clone constructor
      {
        util::ShPtr< assemble::CollectorAASpecified> clone_constr( param_constr.Clone());
        const util::SiPtrList< const biol::AABase> collected_resis( clone_constr->Collect( model));
        BCL_ExampleCheck( collected_resis.GetSize(), 2);
      }

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< assemble::CollectorAASpecified>(), param_constr.GetClassIdentifier());

      // test GetResidueList function
      BCL_ExampleCheck( param_constr.GetResidueList().GetSize(), 3);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test collect function
      {
        const util::SiPtrList< const biol::AABase> collected_resis( param_constr.Collect( model));
        BCL_ExampleCheck( collected_resis.GetSize(), 2);
        BCL_ExampleCheck( ( *collected_resis.Begin())->GetSeqID(), 45);
        BCL_ExampleCheck( ( *collected_resis.Begin())->GetChainID(), 'A');
        BCL_ExampleCheck( ( *++collected_resis.Begin())->GetSeqID(), 100);
        BCL_ExampleCheck( ( *++collected_resis.Begin())->GetChainID(), 'A');
      }

    //////////////////////
    // input and output //
    //////////////////////

      // test Write and read
      {
        WriteBCLObject( param_constr);
        assemble::CollectorAASpecified read_collector;
        ReadBCLObject( read_collector);
        const util::SiPtrList< const biol::AABase> collected_resis( read_collector.Collect( model));
        BCL_ExampleCheck( collected_resis.GetSize(), 2);
      }

    //////////////////////
    // helper functions //
    //////////////////////

      // test ReadAALocators
      {
        // initialize list of locators read from "resi_filename" and created by ReadAALocators
        const storage::List< assemble::LocatorAA> read_locators
        (
          assemble::CollectorAASpecified::ReadAALocators( resi_filename)
        );
        // do checks to make sure function worked
        BCL_ExampleCheck( read_locators.GetSize(), 3);
        BCL_ExampleCheck( read_locators.Begin()->GetAAID(), 45);
        BCL_ExampleCheck( read_locators.Begin()->GetLocatorChain().GetChainID(), 'A');
        BCL_ExampleCheck( ( ++read_locators.Begin())->GetAAID(), 36);
        BCL_ExampleCheck( ( ++read_locators.Begin())->GetLocatorChain().GetChainID(), 'B');
        BCL_ExampleCheck( ( ++++read_locators.Begin())->GetAAID(), 100);
        BCL_ExampleCheck( ( ++++read_locators.Begin())->GetLocatorChain().GetChainID(), 'A');
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleCollectorAASpecified

  const ExampleClass::EnumType ExampleAssembleCollectorAASpecified::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleCollectorAASpecified())
  );
  
} // namespace bcl
