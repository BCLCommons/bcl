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
#include "command/bcl_command_flag_static_and_dynamic.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "density/bcl_density_protein_agreements.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_flag_static_and_dynamic.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommandFlagStaticAndDynamic :
    public ExampleInterface
  {
  public:

    ExampleCommandFlagStaticAndDynamic *Clone() const
    {
      return new ExampleCommandFlagStaticAndDynamic( *this);
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
      command::ParameterCheckEnumerate< density::ProteinAgreements>::s_Instance.IsDefined();
      // parameter for density map and resolution
      util::ShPtr< command::ParameterInterface> sp_mrc_file_param
      (
        new command::Parameter
        (
          "mrc_file",
          "filename of density map to be used",
          command::ParameterCheckFileExistence(),
          ""
        )
      );
      util::ShPtr< command::ParameterInterface> sp_mrc_resolution_param
      (
        new command::Parameter
        (
          "mrc_resolution",
          "resolution of density map to be used",
          command::ParameterCheckRanged< double>( 0.0, 100.0),
          ""
        )
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // a stream for use in our test cases
      std::stringstream error_stream;

      // flag with dynamic parameters and static parameters
      util::ShPtr< command::FlagStaticAndDynamic> density_score_flag
      (
        new command::FlagStaticAndDynamic
        (
          "density_agreement",
          "information required for density agreement score",
          command::Parameter
          (
            "score",
            "choice of agremment score",
            command::ParameterCheckEnumerate< density::ProteinAgreements>()
          ),
          0, density::GetProteinAgreements().GetEnumCount()
        )
      );
      density_score_flag->PushBack( sp_mrc_file_param);
      density_score_flag->PushBack( sp_mrc_resolution_param);

    /////////////////
    // data access //
    /////////////////

      // number of static parameters
      BCL_ExampleCheck( density_score_flag->GetNumberStaticParameters(), 2);

      // number of dynamic parameters
      BCL_ExampleCheck( density_score_flag->GetNumberDynamicParameters(), 0);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // create a sample parameter for use in testing the ReadFromList function
      const storage::Vector< std:: string> density_flag_vector
      (
        storage::Vector< std:: string>::Create
        (
          AddExampleInputPathToFilename( e_Biology, "1ubi_res_6.6voxelsize_2.200Gaussian.mrc"),
          "6.6",
          density::GetProteinAgreements().e_CCC.GetName()
        )
      );

      // make sure that it can read in a set of parameters from a vector with no trouble
      const bool read1_success( density_score_flag->ReadFromList( density_flag_vector, error_stream));
      BCL_ExampleIndirectCheck( read1_success, true, "ReadFromList: " + error_stream.str());

      // number of dynamic parameters
      BCL_ExampleCheck( density_score_flag->GetNumberDynamicParameters(), 1);

      // content of parameters
      BCL_ExampleCheck( density_score_flag->GetParameterList()( 0)->GetValue(), density_flag_vector( 0));
      BCL_ExampleCheck( density_score_flag->GetParameterList()( 1)->GetValue(), density_flag_vector( 1));
      BCL_ExampleCheck( density_score_flag->GetParameterList()( 2)->GetValue(), density_flag_vector( 2));
      BCL_ExampleCheck( density_score_flag->GetDynamicParameterList()( 0)->GetValue(), density_flag_vector( 2));

    //////////////////////
    // input and output //
    //////////////////////

      // check if it can be written and read from file stream
      WriteBCLObject( *density_score_flag);

      // read file back into different parameter
      command::FlagStaticAndDynamic flag_read;
      ReadBCLObject( flag_read);

      // check that written and read flag are the same
      BCL_ExampleIndirectCheck
      (
           density_score_flag->GetName()        == flag_read.GetName()
        && density_score_flag->GetDescription() == flag_read.GetDescription()
        && density_score_flag->GetFlag()        == flag_read.GetFlag()
        && density_score_flag->GetSize()        == flag_read.GetSize(),
        true,
        "writing and reading"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCommandFlagStaticAndDynamic

  const ExampleClass::EnumType ExampleCommandFlagStaticAndDynamic::s_Instance
  (
    GetExamples().AddEnum( ExampleCommandFlagStaticAndDynamic())
  );

} // namespace bcl

