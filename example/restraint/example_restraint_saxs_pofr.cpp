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
#include "restraint/bcl_restraint_sas_pofr.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_saxs_pofr.cpp
  //! @details compute the distance distribution of protein model
  //!
  //! @author putnamdk
  //! @date July 2, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintSasPofR :
    public ExampleInterface
  {
  public:

    ExampleRestraintSasPofR *Clone() const
    {
      return new ExampleRestraintSasPofR( *this);
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
      // Get the Saxs Data and read in the crysol generated data for ubiquitin
      util::Implementation< restraint::SasDensityData> sp_exp_data( new restraint::SasDensityData());
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1ubi_crysol.gnom"));
      sp_exp_data->ReadFromDataFile( read);
      io::File::CloseClearFStream( read);

      // Read in the protein model
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi_bcl.pdb"), biol::GetAAClasses().e_AAComplete)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      restraint::SasPofR density_result;

      density_result.SetExperimentalDensity( sp_exp_data);

    ///////////////
    // operators //
    ///////////////

      // get the calculated curve
      restraint::SasExperimentalAndCalculatedDensity data_sets( density_result( protein_model));

      BCL_ExampleIndirectCheck
      (
        data_sets.GetCalculatedDensity().GetDensityData().GetSize(),
        sp_exp_data->GetDensityData().GetSize(),
        "() operator"
      );

      BCL_ExampleCheckWithinTolerance
      (
        data_sets.GetCalculatedDensity().GetDmax(),
        43.5178,
        0.0001
      );

      BCL_ExampleCheckWithinTolerance
      (
        data_sets.GetExperimentalDensity().GetDmax(),
        45.75,
        0.0001
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintSasPofR

  const ExampleClass::EnumType ExampleRestraintSasPofR::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintSasPofR())
  );

} // namespace bcl
