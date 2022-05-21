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
#include "restraint/bcl_restraint_saxs_data_reduction.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "io/bcl_io_file.h"
#include "restraint/bcl_restraint_sas_density_data.h"
#include "restraint/bcl_restraint_sas_scattering_data.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_saxs_data_reduction.cpp
  //! @details This example tests the reduction of a large SAXS data set to a subset of the set
  //!          The reduction is based on Shannon Sampling and NoisyData Reduction
  //!
  //! @author putnamdk
  //! @date Aug 1, 2013
  //! @remarks status complete
  //! @remarks reviewed by putnamdk on Aug 1, 2013
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintSaxsDataReduction :
    public ExampleInterface
  {
  public:

      ExampleRestraintSaxsDataReduction *Clone() const
    {
      return new ExampleRestraintSaxsDataReduction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class namee
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // Read Protein Model
      // Read in the protein model
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "3ICL.pdb"), biol::GetAAClasses().e_AAComplete)
      );

      BCL_MessageDbg( "Read in Protein model:");

      // Read and store Unprocessed Experimental Data from GNOM
      util::ShPtr< restraint::SasScatteringData> sp_raw_experimental_data( new restraint::SasScatteringData());

      BCL_MessageDbg( "Created Experimental Data object:");

      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "3ICL_SAXS.gnom"));
      sp_raw_experimental_data->ReadFromDataFile( read);
      io::File::CloseClearFStream( read);

      // Create Shared Pointer to Estimated Reduced Data
      util::ShPtr< restraint::SasScatteringData> sp_estimated_reduced_data
      (
        restraint::SaxsDataReduction::SasSignalRecoveryEstimate( sp_raw_experimental_data, 65.46)
      );

      BCL_MessageDbg( "SasSignalRecoveryEstimate finished ");

      // Create new object
      restraint::SasScatteringData estimated_reduced_read;

      // WriteBCLObject( *sp_estimated_reduced_data);

      // Read saved BCL object ( The correct object)
      std::string estimated_reduction_path( AddExampleInputPathToFilename( e_Biology, "saxs_data_estimated_reduction.data"));

      // Read in object in the bcl_object directory
      ReadBCLObjectfromPath( estimated_reduced_read, estimated_reduction_path);

      // Compare original and read in object
      BCL_ExampleCheck( estimated_reduced_read.GetScatteringData(), sp_estimated_reduced_data->GetScatteringData());

      // Create Shared Pointer to Reduced Data Set
      util::ShPtr< restraint::SasScatteringData> sp_reduced_data
      (
        restraint::SaxsDataReduction::SasSignalRecovery( protein_model, 10, sp_raw_experimental_data, 65.46)
      );

      BCL_MessageDbg( "SasSignalRecovery finished ");

      // Create new object
      restraint::SasScatteringData reduced_read;

      //WriteBCLObject( *sp_reduced_data);

      // Read saved BCL object ( The correct object)
      std::string reduction_path( AddExampleInputPathToFilename( e_Biology, "saxs_data_reduction.data"));

      // Read in object in the bcl_object directory
      ReadBCLObjectfromPath( reduced_read, reduction_path);

      // Compare original and read in object
      BCL_ExampleCheck( reduced_read.GetScatteringData(), sp_reduced_data->GetScatteringData());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintSasScatteringData

  const ExampleClass::EnumType ExampleRestraintSaxsDataReduction::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintSaxsDataReduction())
  );

} // namespace bcl
