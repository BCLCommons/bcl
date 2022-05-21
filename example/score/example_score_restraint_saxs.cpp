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
#include "score/bcl_score_restraint_saxs.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "restraint/bcl_restraint_sas_debye.h"
#include "score/bcl_score_sas_type.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_restraint_saxs.cpp
  //! @details Examle file to demonstrate capability to generate RMSD values from Protein Model
  //!
  //! @author putnamdk
  //! @date May 13, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreRestraintSaxs :
    public ExampleInterface
  {
  public:

    ExampleScoreRestraintSaxs *Clone() const
    {
      return new ExampleScoreRestraintSaxs( *this);
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
      // Get the Saxs Data
      // read in crysol generated data for ubiquitin

      util::ShPtr< restraint::SasScatteringData> sp_exp_data( new restraint::SasScatteringData());
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1ubi_saxs_crysol.gnom"));
      sp_exp_data->ReadFromDataFile( read);
      io::File::CloseClearFStream( read);

      // Get the protein model
      // read in the protein model
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"))
      );

      // Tell the protein model how to access SAXS Data
      // store the experimental data (From Crysol)
      restraint::SasDebye saxs_debye( false, false, 1.0, 0.0, true, false, 0.0);
      saxs_debye.SetExperimentalData( sp_exp_data);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // constructor taking parameters
      const std::string scheme( "test_saxs_score");
      score::RestraintSaxs saxs_param_constr
      (
        saxs_debye,
        ( score::SasType( false, score::SasType::e_chi)),
        scheme
      );
      BCL_ExampleIndirectCheck( saxs_param_constr.GetScheme(), scheme, "constructor");

    ///////////////
    // operators //
    ///////////////

      // test () operator
      BCL_ExampleCheckWithinTolerance( saxs_param_constr( protein_model), 1.53551, 0.001);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintSaxs

  const ExampleClass::EnumType ExampleScoreRestraintSaxs::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreRestraintSaxs())
  );

} // namespace bcl
