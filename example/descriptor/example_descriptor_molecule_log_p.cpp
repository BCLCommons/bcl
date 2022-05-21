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
#include "descriptor/bcl_descriptor_molecule_log_p.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_molecule_misc_property.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_linear_least_squares.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_log_p.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 13, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeLogP :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeLogP *Clone() const
    {
      return new ExampleDescriptorMoleculeLogP( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      descriptor::MoleculeLogP logp;

      // copy constructor
      descriptor::MoleculeLogP logp_copy( logp);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( logp.GetAlias(), "LogP");
      BCL_ExampleCheck( logp.GetString(), "LogP");

    ///////////////
    // operators //
    ///////////////

      std::string filename( "small_molecule_descriptors_example_out_high_resolution.sdf");

      // create input stream for reading a smallmolecule ensemble
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, filename));
      // read in ensemble
      chemistry::FragmentEnsemble ensemble( input);
      // close stream
      io::File::CloseClearFStream( input);

      // keep track of the errors and counts of each atom type
      linal::Vector< float> sum_errors_for_atom_type( chemistry::GetAtomTypes().GetEnumCount(), 0.0);
      linal::Vector< size_t> atom_type_count( chemistry::GetAtomTypes().GetEnumCount(), size_t( 0));

      // keep track of every logp calculated with adriana and chemistry so we can run
      // a correlation on them in the end
      storage::Vector< float> adriana_log_ps;
      storage::Vector< float> chemistry_log_ps;

      // make a small molecule property to get the actual xlogp values from the small molecule
      descriptor::MoleculeMiscProperty adriana_xlogp_getter( "XlogP", 1);

      for
      (
        storage::List< chemistry::FragmentComplete>::iterator itr( ensemble.Begin()), itr_end( ensemble.End());
        itr != itr_end;
        ++itr
      )
      {
        // get the (already calculated) logp in chem
        adriana_log_ps.PushBack( adriana_xlogp_getter( *itr)( 0));

        // calculate the log p
        chemistry_log_ps.PushBack( logp( *itr)( 0));
      }

      // record the chi-squared between the values from adriana and the values from chemistry
      const float adriana_to_chemistry_chi_sq
      (
        1.0 - math::LinearLeastSquares::SolutionAndChiSquared
        (
          linal::Matrix< float>( 1, adriana_log_ps.GetSize(), chemistry_log_ps),
          linal::Vector< float>( chemistry_log_ps)
        ).Second()
      );

      // check that chemistry gave similar results to adriana
      if
      (
        !BCL_ExampleIndirectCheck
        (
          adriana_to_chemistry_chi_sq > 0.9,
          true,
          "correlation of adriana and chemistry uncorrected logp (R^2="
          + util::Format()( adriana_to_chemistry_chi_sq) + ")"
        )
      )
      {
        BCL_MessageCrt
        (
          "adriana values: " + util::Format()( adriana_log_ps)
          + "\nChemistry: " + util::Format()( chemistry_log_ps)
        );
      }

      linal::Vector< float> adriana_chemistry_diff_logps( chemistry_log_ps);
      adriana_chemistry_diff_logps -= linal::Vector< float>( adriana_log_ps);

      const float rmsd_mean_molecular_logp
      (
        math::Statistics::Norm
        (
          adriana_chemistry_diff_logps.Begin(),
          adriana_chemistry_diff_logps.End()
        ) / adriana_chemistry_diff_logps.GetSize()
      );

      // check that adriana and chemistry gave similar results
      // the differences are largely because in chemistry we do some extra extrapolation for some of the phosphorus
      // atom types; in adriana most of the phosphorus types were given a 0 value
      // Also, the check for whether a C_TrTrTrPi is branched trigonal or not was implemented incorrectly in chem
      if
      (
        !BCL_ExampleIndirectCheck
        (
          rmsd_mean_molecular_logp < 10.0,
          true,
          "rmsd between logp calculated in adriana vs chemistry was "
          + util::Format()( rmsd_mean_molecular_logp)
        )
      )
      {
        // print out the actual vectors too:
        BCL_MessageDbg
        (
          "Chemistry logps: " + util::Format()( chemistry_log_ps)
          + "adriana logps: " + util::Format()( adriana_log_ps)
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( logp, logp_copy),
        true,
        "logp I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeLogP

  const ExampleClass::EnumType ExampleDescriptorMoleculeLogP::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeLogP())
  );

} // namespace bcl
