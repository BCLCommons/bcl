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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "density/bcl_density_fit_protein_minimizer_mc.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "coord/bcl_coord_move_rotate_random.h"
#include "coord/bcl_coord_move_transform_random.h"
#include "coord/bcl_coord_move_translate_random.h"
#include "density/bcl_density_map.h"
#include "fold/bcl_fold_mutate_protein_model.h"
#include "math/bcl_math_mutate_decision_node.h"
#include "mc/bcl_mc_approximator.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_criterion_unimproved.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> FitProteinMinimizerMC::s_Instance
    (
      GetObjectInstances().AddInstance( new FitProteinMinimizerMC())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FitProteinMinimizerMC::FitProteinMinimizerMC()
    {
    }

    //! @brief Clone function
    //! @return pointer to new FitProteinMinimizerMC
    FitProteinMinimizerMC *FitProteinMinimizerMC::Clone() const
    {
      return new FitProteinMinimizerMC( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FitProteinMinimizerMC::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the resolution of the density map
    //! @param RESOLUTION density map and simulation resolution
    void FitProteinMinimizerMC::SetResolution( const double RESOLUTION)
    {
      m_Resolution = RESOLUTION;
    }

    //! @brief set max translation and rotation
    //! @param MAX_TRANSLATION max translation in any direction for a single iteration
    //! @param MAX_ROTATION max rotation in radians in any direction for a single iteration
    void FitProteinMinimizerMC::SetMaxTranslationAndRotation( const double MAX_TRANSLATION, const double MAX_ROTATION)
    {
      m_MaxTranslation = MAX_TRANSLATION;
      m_MaxRotation    = MAX_ROTATION;
    }

    //! @brief set the max number of iterations for minimization
    //! @param MAX_NUMBER_ITERATIONS maximum number of iterations for minimization
    void FitProteinMinimizerMC::SetMaxIterations( const size_t MAX_NUMBER_ITERATIONS)
    {
      m_MaxIterations = MAX_NUMBER_ITERATIONS;
    }

    //! @brief set protein agreement measure to be used
    //! @param AGREEMENT protein agreement enumerator
    void FitProteinMinimizerMC::SetProteinAgreement( const ProteinAgreement &AGREEMENT)
    {
      m_ProteinAgreement = AGREEMENT;
    }

    //! @brief simulator to use
    //! @param DENSITY_SIMULATOR simulator enumerator
    void FitProteinMinimizerMC::SetSimulator( const Simulator &DENSITY_SIMULATOR)
    {
      m_Simulator = DENSITY_SIMULATOR;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator minimizing the position of a protein model within a given density map
    //! @param PROTEIN_MODEL start position of given protein model
    //! @param DENSITY_MAP the density map to fit the PROTEIN_MODEL into
    //! @return the fitted protein model
    assemble::ProteinModel FitProteinMinimizerMC::operator()( const assemble::ProteinModel &PROTEIN_MODEL, const Map &DENSITY_MAP) const
    {
      // construct a function that calculates the deviation between simulated density from atoms and the given density map
      util::ShPtr< ProteinAgreementInterface> density_agreement
      (
        GetProteinAgreements().CreateProteinAgreement
        (
          m_ProteinAgreement,
          m_Simulator,
          util::ToSiPtr( DENSITY_MAP),
          m_Resolution
        )
      );

      // mutate rotate
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > sp_mutate_rotate
      (
        new fold::MutateProteinModel( coord::MoveRotateRandom( m_MaxRotation))
      );

      // mutate translate
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > sp_mutate_translate
      (
        new fold::MutateProteinModel( coord::MoveTranslateRandom( m_MaxTranslation))
      );

      // mutate function rot trans
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > sp_mutate_rot_trans
      (
        new fold::MutateProteinModel( coord::MoveTransformRandom( m_MaxTranslation, m_MaxRotation))
      );

      math::ObjectProbabilityDistribution< math::MutateInterface< assemble::ProteinModel> > mutates;
      mutates.PushBack( 0.25, *sp_mutate_translate);
      mutates.PushBack( 0.25, *sp_mutate_rotate);
      mutates.PushBack( 0.25, *sp_mutate_rot_trans);

      // mutate function
      math::MutateDecisionNode< assemble::ProteinModel> mutate( mutates);

//      // create printer
//      util::ShPtr< mc::MoviePrinterInterface> sp_movie
//      (
//        new mc::MoviePrinterChimera()
//      );
//
//      // if user wished to print all intermediate minimization steps
//      if( m_WriteMinimizationFlag->GetFlag())
//      {
//        const storage::Set< opti::StepStatusEnum> step_status_set
//        (
//          opti::e_Improved, opti::e_Accepted, opti::e_Skipped, opti::e_Rejected
//        );
//
//        // initialize movie
//        sp_movie->Initialize
//        (
//          m_OutputPrefix->GetFirstParameter()->GetValue(),
//          storage::Vector< std::string>::Create( GetStaticClassName< storage::Table< double> >(), "value"),
//          storage::Vector< std::string>( 1, density_agreement->GetScheme()),
//          720, 480, false
//        );
//
//        // construct movie printer
//        util::ShPtr< assemble::PrinterProteinModelMovie> sp_printer
//        (
//          new assemble::PrinterProteinModelMovie
//          (
//            m_OutputPrefix->GetFirstParameter()->GetValue(),
//            sp_movie,
//            density_agreement,
//            step_status_set,
//            quality::GetSuperimposeMeasures().e_NoSuperimpose,
//            storage::Set< quality::Measure>()
//          )
//        );
//
//        // insert density map as first frame
//        sp_movie->SetStartFrame( m_MrcFilenameParam->GetValue(), true);
//
//        // set the water mark
//        sp_movie->SetWaterMark();
//
//        // set the printer
//        sp_tracker->SetPrinter( sp_printer);
//      }

      // create the temperature control
      util::ShPtr< mc::TemperatureInterface> sp_temperature
      (
        new mc::TemperatureAccepted
        (
          mc::TemperatureAccepted::GetParameterStartFraction()->GetNumericalValue< double>(),
          mc::TemperatureAccepted::GetParameterEndFraction()->GetNumericalValue< double>(),
          m_MaxIterations,
          mc::TemperatureAccepted::GetParameterStartTemperature()->GetNumericalValue< double>(),
          mc::TemperatureAccepted::GetParameterUpdateInterval()->GetNumericalValue< size_t>()
        )
      );

      // create the metropolis
      mc::Metropolis< double> metropolis( sp_temperature, true);

      // create the termination criterion
      opti::CriterionCombine< assemble::ProteinModel, double> termination_criterion;
      termination_criterion.InsertCriteria
      (
        opti::CriterionNumberIterations< assemble::ProteinModel, double>
        (
          m_MaxIterations
        )
      );
      termination_criterion.InsertCriteria
      (
        opti::CriterionUnimproved< assemble::ProteinModel, double>
        (
          m_MaxIterations / s_MaxIterationsStepsInARowFraction
        )
      );

      // create the approximator
      mc::Approximator< assemble::ProteinModel, double> approximator
      (
        *density_agreement,
        mutate,
        metropolis,
        termination_criterion,
        PROTEIN_MODEL
      );

      // approximate
      approximator.Approximate();

      return approximator.GetTracker().GetBest()->First();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FitProteinMinimizerMC::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FitProteinMinimizerMC::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace density
} // namespace bcl
