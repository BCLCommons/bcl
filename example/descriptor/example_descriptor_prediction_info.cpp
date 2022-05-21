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
#include "descriptor/bcl_descriptor_prediction_info.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_iterator.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "sdf/bcl_sdf_molecule_reading_pref.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_prediction_info.cpp
  //!
  //! @author mendenjl
  //! @date Jan 19, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorPredictionInfo :
    public ExampleInterface
  {
  public:

    ExampleDescriptorPredictionInfo *Clone() const
    {
      return new ExampleDescriptorPredictionInfo( *this);
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
      // get the directory containing the models
      const std::string directory( AddExampleInputPathToFilename( e_Descriptor, "models/PredictionInfoTest"));

      // create an implementation to do the prediction statistics
      descriptor::CheminfoProperty prediction_stats
      (
        "PredictionInfo(predictor=File(directory=" + directory + ", prefix=model),statistics(Min,Max,Mean,StandardDeviation))"
      );
      // and one to compute a few metrics: local ppv (likelihood of activity), ppv (cumulative P(activity)), and fraction
      // of actives predicted above this value (TPR)
      descriptor::CheminfoProperty prediction_metrics
      (
        "PredictionInfo(predictor=File(directory=" + directory + ", prefix=model),metrics(LocalPPV, PPV, TPR))"
      );

      // check sizes
      BCL_ExampleCheck( prediction_stats->GetSizeOfFeatures(), 4);
      BCL_ExampleCheck( prediction_metrics->GetSizeOfFeatures(), 3);

      // create input stream for reading a smallmolecule ensemble
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "histidine.sdf"));
      // read in ensemble
      chemistry::FragmentEnsemble ensemble( input, sdf::e_Saturate);
      io::File::CloseClearFStream( input);

      // now retrieve an active molecule. This happens to be the highest predicted active by this naive linear model
      // based only on scalar descriptors
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "1798_active.sdf"));
      ensemble.ReadMoreFromMdl( input, sdf::e_Saturate);
      const chemistry::FragmentComplete &histidine( *ensemble.Begin());
      const chemistry::FragmentComplete &first_active( *++ensemble.Begin());

      // check the descriptions produced
      BCL_ExampleCheckWithinTolerance
      (
        prediction_stats->SumOverObject( histidine),
        linal::MakeVector< float>( 0.00190348, 0.00320547, 0.00243092, 0.000476654),
        0.005
      );
      BCL_ExampleCheckWithinTolerance
      (
        prediction_stats->SumOverObject( first_active),
        linal::MakeVector< float>( -0.024361, -0.00903787, -0.0180861, 0.00629036),
        0.005
      );
      BCL_ExampleCheckWithinTolerance
      (
        prediction_metrics->SumOverObject( histidine),
        linal::MakeVector< float>( 0.00219683, 0.00370257, 0.776471),
        0.005
      );
      BCL_ExampleCheckWithinTolerance
      (
        prediction_metrics->SumOverObject( first_active),
        linal::MakeVector< float>( 0.0, 0.00302489, 1.0),
        0.005
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorPredictionInfo

  const ExampleClass::EnumType ExampleDescriptorPredictionInfo::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorPredictionInfo())
  );

} // namespace bcl
