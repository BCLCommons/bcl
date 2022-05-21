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
#include "nmr/bcl_nmr_residual_dipolar_coupling_least_square_deviation.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_statistics.h"
#include "nmr/bcl_nmr_star_rdc_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_nmr_residual_dipolar_coupling_least_square_deviation.cpp
  //!
  //! @author alexanns, weinerbe
  //! @date Jun 11, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleNmrResidualDipolarCouplingLeastSquareDeviation :
    public ExampleInterface
  {
  public:

    ExampleNmrResidualDipolarCouplingLeastSquareDeviation *Clone() const
    {
      return new ExampleNmrResidualDipolarCouplingLeastSquareDeviation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // read in the protein model
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"))
      );

      // get the aa sequence
      const util::SiPtrVector< const biol::AASequence> aa_seq( protein_model.GetSequences());

      // read in the restraints
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1ubi_rdc_baxbicHN.exp"));
      nmr::StarRDCHandler handler;
      restraint::RDC rdcs_nh( handler.ReadRestraints( read));
      util::ShPtr< restraint::RDC> nh_restraints( rdcs_nh.Clone());
      io::File::CloseClearFStream( read);

      // read in CH rdcs
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1ubi_rdc_baxbicHC.exp"));
      restraint::RDC rdcs_hc( handler.ReadRestraints( read));
      util::ShPtr< restraint::RDC> hc_restraints( rdcs_hc.Clone());
      io::File::CloseClearFStream( read);

      // read in CH and NH rdcs
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1ubi_rdc_baxbicHCHN.exp"));
      // RDCs need to be normalized so set flag
      nmr::StarRDCHandler handler_normalized( "", 0, false, true);
      restraint::RDC rdcs_hchn( handler_normalized.ReadRestraints( read));
      util::ShPtr< restraint::RDC> hchn_restraints( rdcs_hchn.Clone());
      io::File::CloseClearFStream( read);
      BCL_MessageStd( "Read hchn rdcs")

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      nmr::ResidualDipolarCouplingLeastSquareDeviation default_constr;

      // check for N-H RDCs
      nmr::RDCContainer nh_rdc_container( default_constr( nh_restraints->GenerateAssignment( protein_model)));
      storage::Pair< bool, double> rdc_nh_agreement_and_r
      (
        CheckCalculatedRDCS( nh_rdc_container, GetExpectedNHRDCs())
      );
      BCL_MessageStd( "CheckedCalculatedRDCs")
      BCL_Example_Check
      (
        rdc_nh_agreement_and_r.First(),
        "not all calculated N-H rdcs are correct"
      );
      const double expected_r_nh( 0.9707);
      BCL_MessageStd
      (
        "expected N-H Rvalue is " + util::Format()( expected_r_nh) +
        " and calculated Rvalue is " + util::Format()( rdc_nh_agreement_and_r.Second())
      );
      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance( expected_r_nh, rdc_nh_agreement_and_r.Second(), 0.001),
        "expected C-H Rvalue is " + util::Format()( expected_r_nh) + " but calculated Rvalue is "
        + util::Format()( rdc_nh_agreement_and_r.Second())
      );

      // check for C-H RDCs
      nmr::RDCContainer ch_rdc_container( default_constr( hc_restraints->GenerateAssignment( protein_model)));
      storage::Pair< bool, double> rdc_ch_agreement_and_r
      (
        CheckCalculatedRDCS( ch_rdc_container, GetExpectedCHRDCs())
      );
      BCL_Example_Check
      (
        rdc_ch_agreement_and_r.First(),
        "not all calculated C-H rdcs are correct"
      );
      const double expected_r_ch( 0.9623);
      BCL_MessageStd
      (
        "expected C-H Rvalue is " + util::Format()( expected_r_ch) +
        " and calculated Rvalue is " + util::Format()( rdc_ch_agreement_and_r.Second())
      );
      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance( expected_r_ch, rdc_ch_agreement_and_r.Second(), 0.001),
        "expected C-H Rvalue is " + util::Format()( expected_r_ch) + " but calculated Rvalue is "
        + util::Format()( rdc_ch_agreement_and_r.Second())
      );

      // check for simultaneous C-H and N-H RDCs
      nmr::RDCContainer chnh_rdc_container( default_constr( hchn_restraints->GenerateAssignment( protein_model)));
      storage::Pair< bool, double> rdc_chnh_agreement_and_r
      (
        CheckCalculatedRDCS( chnh_rdc_container, GetExpectedCHNHRDCs())
      );
      BCL_Example_Check
      (
        rdc_chnh_agreement_and_r.First(),
        "not all calculated C-H and N-H rdcs are correct"
      );
      const double bcl_expected_r_chnh( 0.9632);
      const double dipocoup_expected_r_chnh( 0.9632); // dipocoup gives 0.9632 correlation
      BCL_MessageStd
      (
        "expected C-H and N-H Rvalue is " + util::Format()( bcl_expected_r_chnh) +
        " and calculated Rvalue is " + util::Format()( rdc_chnh_agreement_and_r.Second())
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( bcl_expected_r_chnh, rdc_chnh_agreement_and_r.Second()),
        "expected C-H and N-H Rvalue is " + util::Format()( bcl_expected_r_chnh) + " but calculated Rvalue is "
        + util::Format()( rdc_chnh_agreement_and_r.Second())
      );
      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance( dipocoup_expected_r_chnh, rdc_chnh_agreement_and_r.Second(), 0.01),
        "expected C-H and N-H Rvalue is " + util::Format()( dipocoup_expected_r_chnh) + " but calculated Rvalue is "
        + util::Format()( rdc_chnh_agreement_and_r.Second())
      );

      return 0;
    }

    double* GetExpectedNHRDCs() const;
    double* GetExpectedCHRDCs() const;
    double* GetExpectedCHNHRDCs() const;

    storage::Pair< bool, double> CheckCalculatedRDCS
    (
      const nmr::RDCContainer &RDC_CONTAINER, double *EXPECTED_RDCS
    ) const;

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleNmrResidualDipolarCouplingLeastSquareDeviation

  const ExampleClass::EnumType ExampleNmrResidualDipolarCouplingLeastSquareDeviation::s_Instance
  (
    GetExamples().AddEnum( ExampleNmrResidualDipolarCouplingLeastSquareDeviation())
  );

  double *ExampleNmrResidualDipolarCouplingLeastSquareDeviation::GetExpectedNHRDCs() const
  {
    static double expected_rdcs[] =
    {
      -10.3497, 7.876, 9.1208, 8.8933, 7.0715, 3.0774, 4.5484, -8.3901, 5.1657, 4.2423, 9.118, 7.9044, 1.2321, 2.7554
        , -7.789, -1.8925, 1.457, 11.0079, 8.6328, 2.7884, 9.1349, 7.6724, 0.5087, 7.4435, 8.724, 5.3857, 3.5668, 10.0064
      , 7.5586, -10.5526, -13.6856, 10.6226, 3.8311, 9.0116, 7.2298, 7.3658, 7.9654, 3.1326, 6.7151, 8.8145, -6.6667, 4.9979
      , 6.8045, -8.8784, -17.4256, -10.3954, 6.4065, 9.7037, 10.5253, 7.3868, 10.8235, 8.315, 4.9683, 4.3825, -11.0986, 9.7671
      , 4.3584, 10.6421, 8.4174, 8.3193, 5.8643, 3.7656, 6.8863, -0.2736, 8.5614, 5.4414, 7.9748, -0.5442
    };

    return expected_rdcs;
  }

  double *ExampleNmrResidualDipolarCouplingLeastSquareDeviation::GetExpectedCHRDCs() const
  {
    static double expected_rdcs[] =
    {
      -1.3161, 2.1178, -2.7377, -1.7931, -3.1508, -2.431, 2.4331, -1.0945, 2.0035, -2.176, -1.5896, -1.2229, -3.0061
      , -1.5455, -1.2826, -0.4093, -2.2525, -1.3157, 3.9794, -3.0773, -1.1311, 1.1789, 0.5998, -2.0435, -2.229, 3.3133
      , -2.9252, -1.0893, -0.5449, 0.3153, -0.011, -1.8165, 1.4142, -2.465, 0.3509, -2.0421, -2.3759, -1.8796, 0.5398
      , 1.2672, 0.9853, -2.5738, 1.3836, 3.5251, -0.7571, 2.326, -1.6855, -0.5571, 0.9504, -2.697, 1.4599, -1.5517
      , -2.4997, 5.0296, -1.2617, -2.3105, -1.4744, -2.8862, -1.122, -2.1658, -2.1959
    };

    return expected_rdcs;
  }

  double *ExampleNmrResidualDipolarCouplingLeastSquareDeviation::GetExpectedCHNHRDCs() const
    {
      static double expected_rdcs[] =
      {
          3.49298,  -7.68623, 9.47597,  6.24047,  10.7319,  8.3894, -9.41371, 2.7504, -7.76524, 6.59805,  5.12467,
          4.34339,  10.2205,  5.25191,  3.36663,  1.88087,  7.26963,  4.1997, -13.4605, 10.8013,  2.82426,  -3.05823,
          -2.11004, 6.83533,  7.17723,  -10.0688, 10.2488,  2.92921,  2.63326,  -1.36661, -0.459813,  5.63416,  -3.85235,
          7.89889,  -0.429536,  6.29469,  7.98315,  5.4168, -1.35944, -3.86653, -4.51058, 8.27538,  -4.06404, -12.3805,
          3.49538,  -7.74669, 5.43534,  2.79744,  -4.0595,  8.99435,  -3.90721, 4.87558,  8.72201,  -16.9948, 4.35763,
          7.4208, 5.43327,  10.095, 3.08871,  7.35164,  6.72923,  -9.34872, 7.99759,  9.04609,  8.88652,  6.95617,
          3.64079,  5.20737,  -7.34415, 5.7572, 4.5199, 9.06537,  8.01189,  1.68342,  3.49226,  -6.67481, -1.50235,
          2.42615,  10.5371,  8.73557,  2.65029,  9.1155, 7.65677,  0.139592, 7.34665,  8.78673,  5.1173, 3.34821,
          9.95271,  7.69192,  -11.2978, -12.6482, 10.4693,  3.4244, 8.49169,  7.26334,  7.46553,  8.20774,  3.86802,
          7.14687,  8.86776,  -6.00909, 5.15369,  6.67901,  -8.56173, -17.3393, -11.3209, 6.28996,  9.15829,  10.1782,
          7.52541,  10.5158,  7.79056,  5.60684,  3.34418,  -11.5929, 9.72166,  3.40375,  10.4868,  8.11511,  8.3533,
          6.41069,  3.24047,  6.44249,  -1.32286, 7.96657,  5.14862,  8.00734,  0.422071
      };
      return expected_rdcs;
    }

  storage::Pair< bool, double> ExampleNmrResidualDipolarCouplingLeastSquareDeviation::CheckCalculatedRDCS
  (
    const nmr::RDCContainer &RDC_CONTAINER, double *EXPECTED_RDCS
  ) const
  {
    const storage::Vector< double> &exp( RDC_CONTAINER.GetExperimentalValues());
    const storage::Vector< double> &thr( RDC_CONTAINER.GetCalculatedlValues());

    bool all_rdcs_correct( true);
    for
    (
      storage::Vector< double>::const_iterator exp_itr( exp.Begin()), exp_itr_end( exp.End()),
        calc_itr( thr.Begin()), calc_itr_end( thr.End());
      exp_itr != exp_itr_end && calc_itr != calc_itr_end; ++exp_itr, ++calc_itr, ++EXPECTED_RDCS
    )
    {
      BCL_MessageDbg
      (
        "Experimental : " + util::Format()( *exp_itr) + " Theoretical : " +
        util::Format()( *calc_itr) + " Expected " + util::Format()( *EXPECTED_RDCS) + " Deviation " +
        util::Format()( *calc_itr - *EXPECTED_RDCS)
      );
      if( !math::EqualWithinTolerance( *EXPECTED_RDCS, *calc_itr, 0.02, 0.05))
      {
        all_rdcs_correct = false;
        BCL_MessageCrt
        (
          " Theoretical : " + util::Format()( *calc_itr) + " != Expected "
          + util::Format()( *EXPECTED_RDCS) + " Deviation : " +
          util::Format()( *calc_itr - *EXPECTED_RDCS)
        );
      }
    }
    const double r( math::Statistics::CorrelationCoefficient( exp.Begin(), exp.End(), thr.Begin(), thr.End()));
    BCL_MessageDbg
    (
      "\nRsquared : \n" +
      util::Format()( math::Statistics::RSquared( exp.Begin(), exp.End(), thr.Begin(), thr.End())) +
      "\nR is " + util::Format()( r)
    );
    return storage::Pair< bool, double>( all_rdcs_correct, r);
  }

} // namespace bcl
