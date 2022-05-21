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
#include "opencl/bcl_opencl_density_fit_protein_minimizer_powell.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_fit_protein_minimizer_powell.h"
#include "density/bcl_density_fit_protein_minimizers.h"
#include "density/bcl_density_map.h"
#include "opti/bcl_opti_approximator_powell.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_number_iterations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  ///////////
  // types //
  ///////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param RESOLUTION the resolution to simulate for
    DensityFitProteinMinimzerPowell::PositionCorrelation::PositionCorrelation()
    {
    }

    //! @brief Clone function
    //! @return pointer to new PositionCorrelation
    DensityFitProteinMinimzerPowell::PositionCorrelation *DensityFitProteinMinimzerPowell::PositionCorrelation::Clone() const
    {
      return new PositionCorrelation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DensityFitProteinMinimzerPowell::PositionCorrelation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the protein model
    //! @param PROTEIN_MODEL
    void DensityFitProteinMinimzerPowell::PositionCorrelation::SetProtein( const assemble::ProteinModel &PROTEIN_MODEL)
    {
      m_Atoms            = m_Simulator.AtomsToDevice( PROTEIN_MODEL.GetAtoms());
      m_NrAtoms          = m_Atoms.GetNumberRows();
    }

    //! @brief set the map
    //! @param DENSITY_MAP
    void DensityFitProteinMinimzerPowell::PositionCorrelation::SetDensityMap( const util::SiPtr< const density::Map> &DENSITY_MAP)
    {
      m_SpMap            = DENSITY_MAP;
      m_Dimensions       = m_Simulator.RoundUpDimensions( m_SpMap->GetDimensions());
      m_DensityMapBuffer = m_Correlation.MapToDevice( *m_SpMap, m_Dimensions);
    }

    //! @brief set the resolution
    //! @param RESOLUTION
    void DensityFitProteinMinimzerPowell::PositionCorrelation::SetResolution( const double RESOLUTION)
    {
      m_Simulator.SetResolution( RESOLUTION);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief is this class compatible with given command queue
    //! @param COMMAND_QUEUE the command queue this object would operate on
    //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
    bool DensityFitProteinMinimzerPowell::PositionCorrelation::IsCompatible( const CommandQueue &COMMAND_QUEUE) const
    {
      return m_Simulator.IsCompatible( COMMAND_QUEUE) &&
             m_Correlation.IsCompatible( COMMAND_QUEUE);
    }

    //! @brief initialize this class
    //! @brief COMMAND_QUEUE queue to use
    //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
    bool DensityFitProteinMinimzerPowell::PositionCorrelation::Initialize( const CommandQueue &COMMAND_QUEUE)
    {
      m_CommandQueue = COMMAND_QUEUE;
      m_MinMax = DataSetMinMax< double>( m_CommandQueue);
      // initialize opencl classes
      return m_Simulator.Initialize( m_CommandQueue) &&
             m_Correlation.Initialize( m_CommandQueue) &&
             m_Transformer.Initialize( m_CommandQueue);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the correlation
    //! @param VECTOR containing 6 elements - three rotations, three translations
    double DensityFitProteinMinimzerPowell::PositionCorrelation::operator()( const linal::Vector< double> &VECTOR) const
    {
      const math::TransformationMatrix3D transformation( VECTOR);

      // transform the coordinates
      const Matrix< double> atoms( m_Transformer( m_Atoms, transformation, m_NrAtoms));

      linal::Vector3D mincoord( m_MinMax.Min( atoms).Begin());
      linal::Vector3D maxcoord( m_MinMax.Max( atoms).Begin());

      // add margin
      mincoord -= 2 * m_SpMap->GetCellWidth();
      maxcoord += 2 * m_SpMap->GetCellWidth();

      // determine index
      const storage::VectorND< 3, int> index
      (
        int( std::floor( mincoord.X() / m_SpMap->GetCellWidth().X())),
        int( std::floor( mincoord.Y() / m_SpMap->GetCellWidth().Y())),
        int( std::floor( mincoord.Z() / m_SpMap->GetCellWidth().Z()))
      );

      // dimensions of grid
      const storage::VectorND< 3, size_t> exact_dimensions
      (
        size_t( std::ceil( maxcoord.X() / m_SpMap->GetCellWidth().X())) - index.First() ,
        size_t( std::ceil( maxcoord.Y() / m_SpMap->GetCellWidth().Y())) - index.Second(),
        size_t( std::ceil( maxcoord.Z() / m_SpMap->GetCellWidth().Z())) - index.Third()
      );
      // dimensions of grid
      const storage::VectorND< 3, size_t> dimensions( m_Simulator.RoundUpDimensions( exact_dimensions));

      // relative index of argument density relative to this
      linal::Vector< int> this_start( 3, int( 0));
      linal::Vector< int> arg_start(  3, int( 0));
      linal::Vector< int> this_end(   3, int( 0));
      linal::Vector< int> arg_end(    3, int( 0));
      // find common start and end
      linal::Vector< int> common_start( 3, int( 0));
      linal::Vector< int> common_end( 3, int( 0));
      linal::Vector< int> extent( 3, int( 0));

      // iterate over dimensions
      for( size_t i( 0); i < 3; ++i)
      {
        this_start( i) = m_SpMap->GetIndex()( i);
        arg_start(  i) =               index( i);
        const int this_end      = this_start( i) + m_Dimensions( i);
        const int arg_end       = arg_start(  i) + dimensions( i);
        const int common_start  = std::max( this_start( i), arg_start( i));
        const int common_end    = std::min( this_end      , arg_end);
        extent( i) = common_end - common_start - 1;

        // no common sub density
        if( extent( i) <= 0)
        {
          return 0.0;
        }

        this_start( i) = common_start - m_SpMap->GetIndex()( i);
        arg_start(  i) = common_start -               index( i);
      }

      // buffer for the density map
      const Vector< double> density( m_Simulator.Simulate( atoms, m_NrAtoms, index, dimensions));

      // calculate cross correlation
      const double ccc
      (
        -m_Correlation.CrossCorrelationCoefficient
        (
          m_DensityMapBuffer,
          density,
          this_start,
          m_Dimensions,
          arg_start,
          dimensions,
          extent,
          0.0
        )
      );

      // end
      return ccc;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DensityFitProteinMinimzerPowell::PositionCorrelation::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DensityFitProteinMinimzerPowell::PositionCorrelation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DensityFitProteinMinimzerPowell::s_Instance
    (
      GetObjectInstances().AddInstance( new DensityFitProteinMinimzerPowell())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DensityFitProteinMinimzerPowell::DensityFitProteinMinimzerPowell() :
      m_CorrelationFunction( new PositionCorrelation())
    {
    }

    //! @brief Clone function
    //! @return pointer to new DensityFitProteinMinimzerPowell
    DensityFitProteinMinimzerPowell *DensityFitProteinMinimzerPowell::Clone() const
    {
      return new DensityFitProteinMinimzerPowell( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DensityFitProteinMinimzerPowell::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the resolution of the density map
    //! @param RESOLUTION density map and simulation resolution
    void DensityFitProteinMinimzerPowell::SetResolution( const double RESOLUTION)
    {
      m_CorrelationFunction->SetResolution( RESOLUTION);
    }

    //! @brief set max translation and rotation
    //! @param MAX_TRANSLATION max translation in any direction for a single iteration
    //! @param MAX_ROTATION max rotation in radians in any direction for a single iteration
    void DensityFitProteinMinimzerPowell::SetMaxTranslationAndRotation( const double MAX_TRANSLATION, const double MAX_ROTATION)
    {
      m_MaxTranslation = MAX_TRANSLATION;
      m_MaxRotation    = MAX_ROTATION;
    }

    //! @brief set the max number of iterations for minimization
    //! @param MAX_NUMBER_ITERATIONS maximum number of iterations for minimization
    void DensityFitProteinMinimzerPowell::SetMaxIterations( const size_t MAX_NUMBER_ITERATIONS)
    {
      m_MaxIterations = MAX_NUMBER_ITERATIONS;
    }

    //! @brief set protein agreement measure to be used
    //! @param AGREEMENT protein agreement enumerator
    void DensityFitProteinMinimzerPowell::SetProteinAgreement( const density::ProteinAgreement &AGREEMENT)
    {
      if( AGREEMENT != density::GetProteinAgreements().e_CCC)
      {
        BCL_MessageStd( "currently, only CCC as ProteinAgreement is supported via Opencl");
      }
    }

    //! @brief simulator to use
    //! @param DENSITY_SIMULATOR simulator enumerator
    void DensityFitProteinMinimzerPowell::SetSimulator( const density::Simulator &DENSITY_SIMULATOR)
    {
      if( DENSITY_SIMULATOR != density::GetSimulators().e_GaussianSphere)
      {
        BCL_MessageStd( "currently, only GaussianSphere as density simulator is supported via Opencl");
      }
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief is this class compatible with given command queue
    //! @param COMMAND_QUEUE the command queue this object would operate on
    //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
    bool DensityFitProteinMinimzerPowell::IsCompatible( const CommandQueue &COMMAND_QUEUE) const
    {
      return m_CorrelationFunction->IsCompatible( COMMAND_QUEUE);
    }

    //! @brief initialize this class
    //! @brief COMMAND_QUEUE queue to use
    //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
    bool DensityFitProteinMinimzerPowell::Initialize( const CommandQueue &COMMAND_QUEUE)
    {
      return m_CorrelationFunction->Initialize( COMMAND_QUEUE);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator minimizing the position of a protein model within a given density map
    //! @param PROTEIN_MODEL start position of given protein model
    //! @param DENSITY_MAP the density map to fit the PROTEIN_MODEL into
    //! @return the fitted protein model
    assemble::ProteinModel DensityFitProteinMinimzerPowell::operator()( const assemble::ProteinModel &PROTEIN_MODEL, const density::Map &DENSITY_MAP) const
    {
      // set members of correlation function
      m_CorrelationFunction->SetProtein( PROTEIN_MODEL);
      m_CorrelationFunction->SetDensityMap( util::ToSiPtr( DENSITY_MAP));

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
        *m_CorrelationFunction, criterion_combine, search_directions, start
      );

      // do the actual approximation
      approximator.Approximate();

      // create the final model
      const util::ShPtr< assemble::ProteinModel> best_model
      (
        density::FitProteinMinimizerPowell::PositionCorrelation::TransformedHardCopy
        (
          PROTEIN_MODEL, approximator.GetTracker().GetBest()->First()
        )
      );

      return *best_model;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DensityFitProteinMinimzerPowell::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DensityFitProteinMinimzerPowell::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DensityFitProteinMinimzerPowellEnumHandler
    //! @brief handler class for adding the operations enum handler
    //! @author woetzen, loweew
    //! @date Mar 14, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DensityFitProteinMinimzerPowellEnumHandler :
      public signal::Slots
    {

    public:

    //////////
    // data //
    //////////

      //! the enum in the density::FitProteinMinimizers
      density::FitProteinMinimizer e_Minimzer;

      //! the only instance of this class
      static const DensityFitProteinMinimzerPowellEnumHandler s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    private:

      //! @brief default constructor
      DensityFitProteinMinimzerPowellEnumHandler();

    public:

    ////////////////
    // operations //
    ////////////////

      //! @brief update the enum with the command queue from the Tools
      //! @param TOOLS the tolls to get the commandqueue from
      void UpdateEnum( Tools &TOOLS);

    }; // template class DensityFitProteinMinimzerPowellEnumHandler

    //! instance of DensityFitProteinMinimzerPowellEnumHandler
    const DensityFitProteinMinimzerPowellEnumHandler DensityFitProteinMinimzerPowellEnumHandler::s_Instance = DensityFitProteinMinimzerPowellEnumHandler();

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DensityFitProteinMinimzerPowellEnumHandler::DensityFitProteinMinimzerPowellEnumHandler() :
      e_Minimzer( density::GetFitProteinMinimizers().AddEnum( "PowellOpencl", util::ShPtr< DensityFitProteinMinimzerPowell>()))
    {
      GetTools().GetQueueUpdateSignal().Connect( this, &DensityFitProteinMinimzerPowellEnumHandler::UpdateEnum);
    }

    //! @brief update the enum with the command queue from the Tools
    //! @param TOOLS the tolls to get the commandqueue from
    void DensityFitProteinMinimzerPowellEnumHandler::UpdateEnum( Tools &TOOLS)
    {
      util::ShPtr< DensityFitProteinMinimzerPowell> sp_minimizer( new DensityFitProteinMinimzerPowell());
      if( !TOOLS.HasCommandQueues())
      {
        *e_Minimzer = util::ShPtr< DensityFitProteinMinimzerPowell>();
        return;
      }

      // try to initialize
      if( sp_minimizer->Initialize( TOOLS.GetFirstCommandQueue()))
      {
        // just update the existing one with the new one
        *e_Minimzer = sp_minimizer;
      }
      else
      {
        BCL_MessageVrb( "unable to initialize enum: OpenCL");
      }
    }

  } // namespace opencl
} // namespace bcl
