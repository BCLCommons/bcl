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
#include "biol/bcl_biol_aa_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_aa_data.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAAData :
    public ExampleInterface
  {
  public:

    ExampleBiolAAData *Clone() const
    {
      return new ExampleBiolAAData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      biol::AAData aa_data_default;

      // construct from all informations: aatype, seqid, pdbid, pdb insertion code, chainid
      biol::AAData aa_data_info( biol::GetAATypes().CYS, 2, -5, 'K', 'D');

      // copy constructor
      biol::AAData aa_data_copy( aa_data_info);

      // clone
      util::ShPtr< util::ObjectInterface> ptr( aa_data_info.Clone());

    /////////////////
    // data access //
    /////////////////

      // check the class identifier
      BCL_Example_Check
      (
        GetStaticClassName( aa_data_default) == "bcl::biol::AAData",
        "unexpected static class name: " + GetStaticClassName( aa_data_default) + " should be: bcl::biol::AAData"
      );

      // class identifier
      BCL_Example_Check
      (
        ptr->GetClassIdentifier() == GetStaticClassName< biol::AAData>(),
        "unexpected class identifier class name: " + ptr->GetClassIdentifier() + " should be: " + GetStaticClassName< biol::AAData>()
      );

      // aatype
      BCL_MessageStd( "type of this AAData: " + aa_data_info.GetType().GetName());
      BCL_Example_Check
      (
        aa_data_info.GetType() == biol::GetAATypes().CYS,
        "constructed aa_data should be CYS"
      );

      // seqid
      BCL_MessageStd( "seqid of this AAData: " + util::Format()( aa_data_info.GetSeqID()));
      BCL_Example_Check
      (
        aa_data_info.GetSeqID() == 2,
        "constructed aa_data should have seqid 2"
      );

      // pdbid
      BCL_MessageStd( "pdbid of this AAData: " + util::Format()( aa_data_info.GetPdbID()));
      BCL_Example_Check
      (
        aa_data_info.GetPdbID() == -5,
        "constructed aa_data should have pdbid -5"
      );

      // pdbicode
      BCL_MessageStd( "pdbicode of this AAData: " + util::Format()( aa_data_info.GetPdbICode()));
      BCL_Example_Check
      (
        aa_data_info.GetPdbICode() == 'K',
        "constructed aa_data should have pdb insertion code K"
      );

      // chain id
      BCL_MessageStd( "chain id of this AAData: " + util::Format()( aa_data_info.GetChainID()));
      BCL_Example_Check
      (
        aa_data_info.GetChainID() == 'D',
        "constructed aa_data should have chainid D"
      );

      // blast profile
      BCL_MessageStd( "blast profile of this AAData: " + util::Format()( aa_data_info.GetBlastProfile()));
      BCL_Example_Check
      (
        !aa_data_info.GetBlastProfile().IsDefined(),
        "constructed aa_data should have undefined blast profile"
      );

      // ss predications as map
      BCL_Example_Check
      (
        aa_data_info.GetSSPredictions().IsEmpty(),
        "constructed aa_data should have empty ss predictions map"
      );

      // particular ss prediction
      BCL_Example_Check
      (
        !aa_data_info.GetSSPrediction( sspred::GetMethods().e_JUFO).IsDefined(),
        "constructed aa_data should not have JUFO ss prediction"
      );

      // check members of default AAData
      BCL_Example_Check
      (
        !aa_data_default.GetType().IsDefined() &&
        aa_data_default.GetSeqID() == 1 &&
        aa_data_default.GetPdbID() == 1 &&
        aa_data_default.GetPdbICode() == ' ' &&
        aa_data_default.GetChainID() == 'A' &&
        !aa_data_default.GetBlastProfile().IsDefined() &&
        aa_data_default.GetSSPredictions().IsEmpty(),
        "constructed aa_data_default is incorrect: " + util::Format()( aa_data_default)
      );

    ///////////////
    // operators //
    ///////////////

      // smaller than
      BCL_Example_Check
      (
        aa_data_default < aa_data_info,
        "default constructed aadata should come before constructed aa data"
      );
      BCL_Example_Check
      (
        aa_data_default < aa_data_info,
        "default constructed aadata should come before constructed aa data"
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities for biol::AAData");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( aa_data_info);
      BCL_MessageVrb( "read object");
      biol::AAData aa_data_read;
      ReadBCLObject( aa_data_read);

      // compare the objects
      BCL_Example_Check
      (
        aa_data_info.GetType() == aa_data_read.GetType() &&
        aa_data_info.GetSeqID() == aa_data_read.GetSeqID() &&
        aa_data_info.GetPdbID() == aa_data_read.GetPdbID() &&
        aa_data_info.GetPdbICode() == aa_data_read.GetPdbICode() &&
        aa_data_info.GetChainID() == aa_data_read.GetChainID() &&
        aa_data_info.GetBlastProfile().IsDefined() == aa_data_read.GetBlastProfile().IsDefined() &&
        aa_data_info.GetSSPredictions().IsEmpty() == aa_data_read.GetSSPredictions().IsEmpty(),
        "read aadata is different from written"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolAAData

  const ExampleClass::EnumType ExampleBiolAAData::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAAData())
  );

} // namespace bcl

