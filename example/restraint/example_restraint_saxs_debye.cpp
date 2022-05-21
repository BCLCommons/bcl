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
#include "restraint/bcl_restraint_sas_debye.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_saxs_debye.cpp
  //! @details compare the calculated intensity with experimental intensity from CRYSOL for ubiquitin
  //!
  //! @author putnamdk
  //! @date May 9, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintSasDebye :
    public ExampleInterface
  {
  public:

    ExampleRestraintSasDebye *Clone() const
    {
      return new ExampleRestraintSasDebye( *this);
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
      util::ShPtr< restraint::SasScatteringData> sp_exp_data( new restraint::SasScatteringData());
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1ubi_saxs_crysol.gnom"));
      sp_exp_data->ReadFromDataFile( read);
      io::File::CloseClearFStream( read);

      // Read in the protein model
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"))
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      restraint::SasDebye result;
      result.SetExperimentalData( sp_exp_data);

    ///////////////
    // operators //
    ///////////////

      // get the calculated curve
      restraint::SasExperimentalAndCalculatedData data_sets( result( protein_model));
      BCL_ExampleIndirectCheck
      (
        data_sets.GetCalculatedData().GetScatteringData().GetSize(),
        sp_exp_data->GetScatteringData().GetSize(),
        "() operator"
      );

      // iterate over experimental and calculated data to get values
      data_sets.WriteToGnuplot( util::GetLogger());

      // test Sans Implementation
      restraint::SasDebye sans_result( false, false, 1.0, 0.0, false, true, 1.0);
      sans_result.SetExperimentalData( sp_exp_data);

      // get the calculated curve
      restraint::SasExperimentalAndCalculatedData sans_data_sets( sans_result( protein_model));

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintSasDebye

  const ExampleClass::EnumType ExampleRestraintSasDebye::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintSasDebye())
  );

} // namespace bcl
