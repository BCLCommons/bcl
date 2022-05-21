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
#include "assemble/bcl_assemble_collector_aa_type.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_collector_aa_type.cpp
  //!
  //! @author alexanns
  //! @date Mar 3, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleCollectorAAType :
    public ExampleInterface
  {
  public:

    ExampleAssembleCollectorAAType *Clone() const
    {
      return new ExampleAssembleCollectorAAType( *this);
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

      // create proteins
      const util::ShPtr< assemble::ProteinModel> ubi
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb")).Clone()
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      assemble::CollectorAAType def_constr;
      util::SiPtrList< const biol::AABase> collected_aas( def_constr.Collect( ubi->GetAminoAcids()));
      BCL_MessageDbg
      (
        "default constructor collected " + util::Format()( collected_aas.GetSize()) +
        " residues. They are " + util::Format()( collected_aas)
      );
      BCL_ExampleCheck( collected_aas.GetSize() == 6, true);

      // constructor taking parameters
      assemble::CollectorAAType param_const( storage::Set< biol::AAType>( biol::GetAATypes().GLY));
      BCL_ExampleCheck( param_const.Collect( ubi->GetAminoAcids()).GetSize(), 6);

      // clone constructor
      util::ShPtr< assemble::CollectorAAType> clone_constr( param_const.Clone());
      BCL_ExampleCheck( param_const.Collect( ubi->GetAminoAcids()).GetSize(), 6);

    /////////////////
    // data access //
    /////////////////

      // GetStaticClassName
      const std::string correct_static_class_name( "bcl::assemble::CollectorAAType");
      BCL_ExampleCheck( GetStaticClassName< assemble::CollectorAAType>(), correct_static_class_name);

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< assemble::CollectorAAType>(), clone_constr->GetClassIdentifier());

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // collect
      {
        assemble::CollectorAAType test( storage::Set< biol::AAType>::Create( biol::GetAATypes().GLY, biol::GetAATypes().LYS));
        BCL_ExampleCheck( test.Collect( ubi->GetAminoAcids()).GetSize(), 13);
      }

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      // write the object
      WriteBCLObject( param_const);

      // read the object back in
      assemble::CollectorAAType read;
      ReadBCLObject( read);
      BCL_ExampleCheck( param_const.Collect( ubi->GetAminoAcids()).GetSize(), 6);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleCollectorAAType

  const ExampleClass::EnumType ExampleAssembleCollectorAAType::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleCollectorAAType())
  );
  
} // namespace bcl
