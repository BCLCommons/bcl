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
#include "density/bcl_density_fit_protein_minimizer_powell.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "density/bcl_density_map.h"
#include "opti/bcl_opti_approximator_powell.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_number_iterations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

  ///////////
  // types //
  ///////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FitProteinMinimizerPowell::PositionCorrelation::PositionCorrelation
    (
      const util::ShPtr< ProteinAgreementInterface> &AGREEMENT,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) :
      m_Agreement( AGREEMENT),
      m_Model( util::ToSiPtr( PROTEIN_MODEL))
    {
    }

    //! @brief Clone function
    //! @return pointer to new PositionCorrelation
    FitProteinMinimizerPowell::PositionCorrelation *FitProteinMinimizerPowell::PositionCorrelation::Clone() const
    {
      return new PositionCorrelation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FitProteinMinimizerPowell::PositionCorrelation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the correlation
    //! @param VECTOR containing 6 elements - three rotations, three translations
    double FitProteinMinimizerPowell::PositionCorrelation::operator()( const linal::Vector< double> &VECTOR) const
    {
      // return the agreement
      return m_Agreement->operator ()( *TransformedHardCopy( *m_Model, VECTOR));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FitProteinMinimizerPowell::PositionCorrelation::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FitProteinMinimizerPowell::PositionCorrelation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    util::ShPtr< assemble::ProteinModel> FitProteinMinimizerPowell::PositionCorrelation::TransformedHardCopy
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const linal::Vector< double> &VECTOR
    )
    {
      // copy the protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.HardCopy());

      // create transformation
      const math::TransformationMatrix3D transformation( VECTOR);

      // apply transformation
      new_model->Transform( transformation);

      // end
      return new_model;
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> FitProteinMinimizerPowell::s_Instance
    (
      GetObjectInstances().AddInstance( new FitProteinMinimizerPowell())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FitProteinMinimizerPowell::FitProteinMinimizerPowell()
    {
    }

    //! @brief Clone function
    //! @return pointer to new FitProteinMinimizerPowell
    FitProteinMinimizerPowell *FitProteinMinimizerPowell::Clone() const
    {
      return new FitProteinMinimizerPowell( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FitProteinMinimizerPowell::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the resolution of the density map
    //! @param RESOLUTION density map and simulation resolution
    void FitProteinMinimizerPowell::SetResolution( const double RESOLUTION)
    {
      m_Resolution = RESOLUTION;
    }

    //! @brief set max translation and rotation
    //! @param MAX_TRANSLATION max translation in any direction for a single iteration
    //! @param MAX_ROTATION max rotation in radians in any direction for a single iteration
    void FitProteinMinimizerPowell::SetMaxTranslationAndRotation( const double MAX_TRANSLATION, const double MAX_ROTATION)
    {
      m_MaxTranslation = MAX_TRANSLATION;
      m_MaxRotation    = MAX_ROTATION;
    }

    //! @brief set the max number of iterations for minimization
    //! @param MAX_NUMBER_ITERATIONS maximum number of iterations for minimization
    void FitProteinMinimizerPowell::SetMaxIterations( const size_t MAX_NUMBER_ITERATIONS)
    {
      m_MaxIterations = MAX_NUMBER_ITERATIONS;
    }

    //! @brief set protein agreement measure to be used
    //! @param AGREEMENT protein agreement enumerator
    void FitProteinMinimizerPowell::SetProteinAgreement( const ProteinAgreement &AGREEMENT)
    {
      m_ProteinAgreement = AGREEMENT;
    }

    //! @brief simulator to use
    //! @param DENSITY_SIMULATOR simulator enumerator
    void FitProteinMinimizerPowell::SetSimulator( const Simulator &DENSITY_SIMULATOR)
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
    assemble::ProteinModel FitProteinMinimizerPowell::operator()( const assemble::ProteinModel &PROTEIN_MODEL, const Map &DENSITY_MAP) const
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

      // create function opencl
      util::ShPtr< math::FunctionInterfaceSerializable< linal::Vector< double>, double> > sp_objective
      (
        new PositionCorrelation
        (
          density_agreement,
          PROTEIN_MODEL
        )
      );

      storage::Vector< linal::Vector< double> > search_directions;
      linal::Vector< double> direction( 6, 0.0);
      const linal::Vector< double> start( direction);
      direction( 0) = m_MaxRotation;
      search_directions.PushBack( direction);
      direction( 0) = 0.0;
      direction( 1) = m_MaxRotation;
      search_directions.PushBack( direction);
      direction( 1) = 0.0;
      direction( 2) = m_MaxRotation * 0.5;
      search_directions.PushBack( direction);
      direction( 2) = 0.0;
      direction( 3) = m_MaxTranslation;
      search_directions.PushBack( direction);
      direction( 3) = 0.0;
      direction( 4) = m_MaxTranslation;
      search_directions.PushBack( direction);
      direction( 4) = 0.0;
      direction( 5) = m_MaxTranslation;
      search_directions.PushBack( direction);

      // create termination criteria for the approximation
      opti::CriterionCombine< linal::Vector< double>, double> criterion_combine;
      criterion_combine.InsertCriteria
      (
        opti::CriterionNumberIterations< linal::Vector< double>, double>( m_MaxIterations / 10)
      );
      criterion_combine.InsertCriteria
      (
        opti::CriterionConvergenceResult< linal::Vector< double>, double>( 1, 0.001)
      );

      // create powell approximator from its members
      opti::ApproximatorPowell< linal::Vector< double>, double> approximator
      (
        *sp_objective, criterion_combine, search_directions, start
      );

      // do the actual approximation
      approximator.Approximate();

      // create the final model
      const util::ShPtr< assemble::ProteinModel> best_model
      (
        PositionCorrelation::TransformedHardCopy( PROTEIN_MODEL, approximator.GetTracker().GetBest()->First())
      );

      return *best_model;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FitProteinMinimizerPowell::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FitProteinMinimizerPowell::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace density
} // namespace bcl
