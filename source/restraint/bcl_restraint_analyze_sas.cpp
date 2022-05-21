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
#include "restraint/bcl_restraint_analyze_sas.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "biol/bcl_biol_sasa_data.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram.h"
#include "restraint/bcl_restraint_sas_analysis.h"
#include "restraint/bcl_restraint_sas_debye_interface.h"
#include "restraint/bcl_restraint_sas_optimization.h"
#include "restraint/bcl_restraint_sas_transformation.h"
#include "restraint/bcl_restraint_saxs_opti_result.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeSas::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeSas())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeSas::AnalyzeSas() :
        m_OutFilePostFix( ".AnalyzeSas"),
        m_ExperimentalDataFileName(),
        m_SasaDataFileName(),
        m_ComputedDataFileName(),
        m_ExcludedVolumeParameter( 1.0),
        m_HydrationShellThickness( 0.0),
        m_MaximumDimension( util::GetUndefined< double>()),
        m_OptimizeHydrationShellParameters( false),
        m_C1Min( SasOptimization::default_C1Min),
        m_C1Max( SasOptimization::default_C1Max),
        m_C2Min( SasOptimization::default_C2Min),
        m_C2Max( SasOptimization::default_C2Max),
        m_C1StepSize( SasOptimization::default_C1StepSize),
        m_C2StepSize( SasOptimization::default_C2StepSize),
        m_ScoreType( "chi"),
        m_Cpu( false),
        m_ShouldApproximateSideChains( true),
        m_ShouldApproximateLoops( false),
        m_ShouldUseSans( false),
        m_DeuteriumExchangeParameter( 0.0)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeAtomDistanceHeatmap
    AnalyzeSas *AnalyzeSas::Clone() const
    {
      return new AnalyzeSas( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeSas::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeSas::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief gives the experimental data file name
    //! @return the filename string
    const std::string &AnalyzeSas::GetExperimentalDataFileName() const
    {
      return m_ExperimentalDataFileName;
    }

    //! @brief gives the sasa data file name
    //! @return the filename string
    const std::string &AnalyzeSas::GetSasaDataFileName() const
    {
      return m_SasaDataFileName;
    }

    //! @brief gives the computed data file name
    //! @return the filename string
    const std::string &AnalyzeSas::GetComputedDataFileName() const
    {
      return m_ComputedDataFileName;
    }

    //! @brief gives the excluded volume parameter
    //! @return the c1 parameter
    const float &AnalyzeSas::GetExcludedVolumeParameter() const
    {
      return m_ExcludedVolumeParameter;
    }

    //! @brief gives the hydration shell thickness parameter
    //! @return the c2 parameter
    const float &AnalyzeSas::GetHydrationShellParameter() const
    {
      return m_HydrationShellThickness;
    }

    //! @brief gives the maximum dimension parameter
    //! @return dmax
    const float &AnalyzeSas::GetMaximumDimension() const
    {
      return m_MaximumDimension;
    }

    //! @brief gives the Minimum value for C1 on the search grid
    //! @return C1Min
    const float &AnalyzeSas::GetC1Min() const
    {
      return m_C1Min;
    }

    //! @brief gives the Maximum value for C1 on the search grid
    //! @return C1Max
    const float &AnalyzeSas::GetC1Max() const
    {
      return m_C1Max;
    }

    //! @brief gives the Minimum value for C2 on the search grid
    //! @return C2Min
    const float &AnalyzeSas::GetC2Min() const
    {
      return m_C2Min;
    }

    //! @brief gives the Maximum value for C2 on the search grid
    //! @return C2Max
    const float &AnalyzeSas::GetC2Max() const
    {
      return m_C1Max;
    }

    //! @brief gives the StepSize for C1 on the search grid
    //! @return C1StepSize
    const float &AnalyzeSas::GetC1StepSize() const
    {
      return m_C1StepSize;
    }

    //! @brief gives the StepSize for C2 on the search grid
    //! @return C2StepSize
    const float &AnalyzeSas::GetC2StepSize() const
    {
      return m_C2StepSize;
    }

    const float &AnalyzeSas::GetDeuteriumExchangeParameter() const
    {
      return m_DeuteriumExchangeParameter;
    }

    //! @brief gives name of the score to be used to compare two saxs profiles
    //! @return Score type
    const std::string &AnalyzeSas::GetScoreTypeName() const
    {
      return m_ScoreType;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeSas::GetAlias() const
    {
      static const std::string s_Name( "AnalyzeSas");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzeSas::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {

      // This will need to find a new place to reside.  This provides the ability to read a raw
      // ExperimentalandCalculated Data file and perform any transformation desired without re-computing the
      // sas profiles

      BCL_Assert( !m_ExperimentalDataFileName.empty(), "experimental data file name is not defined");

      // Read and store the Experimental Data.  The SAS curve generated from the protein model will be fit to this data
      util::ShPtr< SasScatteringData> sp_experimental_data( new SasScatteringData());

      io::IFStream read;
      io::File::MustOpenIFStream( read, m_ExperimentalDataFileName);
      sp_experimental_data->ReadFromDataFile( read);
      io::File::CloseClearFStream( read);

      double max_intensity( SasAnalysis::MaxIntensity( *sp_experimental_data));

      if( !m_ComputedDataFileName.empty())
      {
        // Read in ExperimentalandCalculated Data file
        util::ShPtr< SasExperimentalAndCalculatedData> sp_SasExperimentalAndCalculatedData
        (
          new SasExperimentalAndCalculatedData()
        );

        io::IFStream read;
        io::File::MustOpenIFStream( read, m_ComputedDataFileName);
        sp_SasExperimentalAndCalculatedData->ReadFromDataFile( read);
        io::File::CloseClearFStream( read);

        // Transform the data
        SasExperimentalAndCalculatedData transformed_data( m_Transform( *sp_SasExperimentalAndCalculatedData));

        // Score the transformed data
        double result
        (
          score::SasType( m_Transform.GetUseErrors(), m_ScoreType)( transformed_data)
        );

        const storage::Vector< SasTransformation::TransformationTypeEnum> transforms_performed
        (
          m_Transform.GetTransformationTypes()
        );

        std::string transforms;

        for
        (
          storage::Vector< SasTransformation::TransformationTypeEnum>::const_iterator
           transform_itr( transforms_performed.Begin()),
           transform_itr_end( transforms_performed.End());
           transform_itr != transform_itr_end;
          ++transform_itr
        )
        {
          transforms += " " + transform_itr->GetString();
        }

        std::string summary
        (
          "\n ScoreFunction: " + util::Format()( m_ScoreType) + "\n"
          " Max Intensity: " + util::Format()( max_intensity) + "\n"
          " Use errors: " + util::Format()( m_Transform.GetUseErrors()) + "\n"
          " Transforms: " + util::Format()( transforms) + "\n"
          " Score: " + util::Format()( result) + "\n"
          " Hardware: " + util::Format()( m_Cpu) + "\n"
        );

        BCL_MessageStd( summary);

        return util::Format()( result);
      }

      double result( 0.0);

      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

      // iterate over all models in the pdbs.ls file
      for
      (
        assemble::ProteinEnsemble::const_iterator
        model_itr( ENSEMBLE.Begin()), model_itr_end( ENSEMBLE.End());
        model_itr != model_itr_end;
        ++model_itr
      )
      {
        util::ShPtr< assemble::ProteinModel> sp_protein_model( *model_itr);
        assemble::ProteinModel protein_model( *sp_protein_model);

        // Get a shared pointer to the protein model data
        util::ShPtr< assemble::ProteinModelData> sp_protein_model_data( protein_model.GetProteinModelData());

        // If Sasa File Provided insert SASA data into Protein Model
        if( !m_SasaDataFileName.empty())
        {
          // Sasa file must be preprocessed from PDBConvert and MSMS
          util::ShPtr< biol::SasaData> sp_sasa_data( new biol::SasaData());

          io::IFStream read_sasa;
          io::File::MustOpenIFStream( read_sasa, m_SasaDataFileName);
          sp_sasa_data->ReadFromDataFile( read_sasa);
          io::File::CloseClearFStream( read_sasa);

          // Insert the Sasa Data into the Protein Model Data
          sp_protein_model_data->Insert( assemble::ProteinModelData::e_Sasa, sp_sasa_data);
        }

        //Insert the Protein Model Data into the protein model
        protein_model.SetProteinModelData( sp_protein_model_data);

        // Create Local Variable to hold C1 and C2 values.  If The values are optimized, update the variable
        float c1_value( m_ExcludedVolumeParameter), c2_value( m_HydrationShellThickness);

        if( m_OptimizeHydrationShellParameters)
        {
          BCL_Assert( !m_SasaDataFileName.empty(), "sasa data file must be defined to optimize parameters");

          // Create default object
          SasOptimization optimization_object;

          if( !m_DefaultSearchGrid)
          {
            // setup grid optmization with the score type and the error flag
            optimization_object.SetC1Min( m_C1Min);
            optimization_object.SetC1Max( m_C1Max);
            optimization_object.SetC2Min( m_C2Min);
            optimization_object.SetC2Max( m_C2Max);
            optimization_object.SetC1StepSize( m_C1StepSize);
            optimization_object.SetC2StepSize( m_C2StepSize);
          }
          optimization_object.SetScoreFunction( m_ScoreType);
          optimization_object.SetErrorFlag( m_Transform.GetUseErrors());
          optimization_object.SetTransformationTypes( m_Transform);
          optimization_object.SetExperimentalData( sp_experimental_data);
          optimization_object.SetApproximateSideChainsFlag( m_ShouldApproximateSideChains);
          optimization_object.SetHardwareType( m_Cpu);
          optimization_object.SetSasType( m_ShouldUseSans);
          optimization_object.SetDeuteriumExchangeParameter( m_DeuteriumExchangeParameter);

          // Optimize c1 and c2 values. The experimental and sasa data were previously inserted into the protein_model
          SaxsOptiResult optimal_parameters( optimization_object( protein_model));

          // Update c1 and c2 values with optimized parameters
          c1_value = optimal_parameters.GetC1();
          c2_value = optimal_parameters.GetC2();
        }

        // Setup Commandline Strings for either the opencl or non-opencl version of the code
        util::Implementation< SasDebyeInterface> saxs
        (
          SasAnalysis::SetDebyeImplementation
          (
            m_ShouldApproximateLoops,
            m_ShouldApproximateSideChains,
            c1_value,
            c2_value,
            m_Cpu,
            m_ShouldUseSans,
            m_DeuteriumExchangeParameter
          )
        );

        // Set the Experimental Data
        saxs->SetExperimentalData( sp_experimental_data);

        //Create container to hold experimental and computed SAXS profiles and store results in container
        util::ShPtr< SasExperimentalAndCalculatedData> sp_experimental_and_calculated_data
        (
          new SasExperimentalAndCalculatedData( saxs->operator()( protein_model))
        );

        // Transform the data
        SasExperimentalAndCalculatedData transformed_data( m_Transform( *sp_experimental_and_calculated_data));

        // Score the transformed data
        result = score::SasType( m_Transform.GetUseErrors(), m_ScoreType)( transformed_data);

        const storage::Vector< SasTransformation::TransformationTypeEnum> transforms_performed
        (
          m_Transform.GetTransformationTypes()
        );

        std::string transforms;

        for
        (
          storage::Vector< SasTransformation::TransformationTypeEnum>::const_iterator
            transform_itr( transforms_performed.Begin()),
            transform_itr_end( transforms_performed.End());
            transform_itr != transform_itr_end;
          ++transform_itr
        )
        {
          transforms += " " + transform_itr->GetString();
        }

        // get the name of the pdb the structure came from
        util::ShPtr< util::Wrapper< std::string> > pdb_filename
        (
           ( *model_itr)->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
        );

        std::string pdb( "dummy.pdb");

        if( pdb_filename.IsDefined())
        {
          pdb = std::string( pdb_filename->GetData());
        }

        std::string summary
        (
          "\n ScoreFunction: " + util::Format()( m_ScoreType) + "\n"
          " Max_Intensity: " + util::Format()( max_intensity) + "\n"
          " Use errors: " + util::Format()( m_Transform.GetUseErrors()) + "\n"
          " C1: " + util::Format()( c1_value) + "\n"
          " C2: " + util::Format()( c2_value) + "\n"
          " Transforms: " + util::Format()( transforms) + "\n"
          " Score: " + util::Format()( result) + " " + util::Format()( pdb) + "\n"
          " Hardware: " + util::Format()( m_Cpu) + "\n"
        );

        BCL_MessageStd( summary);
      }

      return util::Format()( result);

      // return the stringSasAnalysis
      //return results;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief hook that derived classes can override, called after TryRead successfully reads a serializer containing only initializer info (no data variables)
    //! @param SERIALIZER the serializer object with initialization information
    //! @param ERR_STREAM stream to write out errors to
    //! @return true, unless there were new errors
    bool AnalyzeSas::ReadInitializerSuccessHook( const util::ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
    {
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeSas::GetSerializer() const
    {
      io::Serializer parameters( m_Transform.GetSerializer());
      parameters.SetClassDescription
      (
        "Computes a SAXS profile from PDB or BCL model and compares the curve with Experimental SAXS Data."
        "hydration layer controlled by excluded solvent c1, and border thickness c2 variables."
        "The Solvent assessible surface area (SASA) file is read in from MSMS"
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AnalyzeSas"
      );

      parameters.AddInitializer
      (
        "experimental_profile",
        "the experimental SAXS data",
        io::Serialization::GetAgentInputFilename( &m_ExperimentalDataFileName)
      );

      parameters.AddOptionalInitializer
      (
        "sasa_profile",
        "the solvent accessible surface area of the protein model",
        io::Serialization::GetAgentInputFilename( &m_SasaDataFileName)
      );

      parameters.AddOptionalInitializer
      (
        "computed_profile",
        "the computed SAXS data",
        io::Serialization::GetAgentInputFilename( &m_ComputedDataFileName)
      );

      parameters.AddInitializer
      (
        "c1",
        "the excluded volume tuning parameter",
        io::Serialization::GetAgent( &m_ExcludedVolumeParameter)
      );

      parameters.AddInitializer
      (
        "c2",
        "the hydration shell thickness",
        io::Serialization::GetAgent( &m_HydrationShellThickness)
      );

      parameters.AddOptionalInitializer
      (
        "c1_min",
        "the minimum excluded volume parameter value on the grid search",
        io::Serialization::GetAgent( &m_C1Min)
      );

      parameters.AddOptionalInitializer
      (
        "c1_max",
        "the maximum excluded volume parameter value on the grid search",
        io::Serialization::GetAgent( &m_C1Max)
      );

      parameters.AddOptionalInitializer
      (
        "c2_min",
        "the minimum hydration thickness parameter value on the grid search",
        io::Serialization::GetAgent( &m_C2Min)
      );

      parameters.AddOptionalInitializer
      (
        "c2_max",
        "the maximum hydration thickness parameter value on the grid search",
        io::Serialization::GetAgent( &m_C2Max)
      );

      parameters.AddOptionalInitializer
      (
        "c1_stepsize",
        "the excluded volume parameter stepsize on the grid search",
        io::Serialization::GetAgent( &m_C1StepSize)
      );

      parameters.AddOptionalInitializer
      (
        "c2_stepsize",
        "the hydration shell parameter stepsize on the grid search",
        io::Serialization::GetAgent( &m_C2StepSize)
      );

      parameters.AddOptionalInitializer
      (
        "dmax",
        "the maximum dimension of the particle from Gnome",
        io::Serialization::GetAgent( &m_MaximumDimension)
      );

      parameters.AddOptionalInitializer
      (
        "optimize_hydration_parameters",
        "1 if the hydration parameters should be optimized - 0 otherwise",
        io::Serialization::GetAgent( &m_OptimizeHydrationShellParameters)
      );

      parameters.AddInitializer
      (
        "default_search_grid",
        "1 to use default parameters- 0 otherwise",
        io::Serialization::GetAgent( &m_DefaultSearchGrid)
      );

      parameters.AddInitializer
      (
        "scoring_function",
        "The scoring function to use for analysis",
        io::Serialization::GetAgent( &m_ScoreType)
      );

      parameters.AddOptionalInitializer
       (
         "cpu",
         "true to force cpu implementation of debye formula - false otherwise",
         io::Serialization::GetAgent( &m_Cpu)
       );

      parameters.AddInitializer
       (
         "approximate_side_chains",
         "true sum form factors to cb atom - false otherwise",
         io::Serialization::GetAgent( &m_ShouldApproximateSideChains)
       );

      parameters.AddInitializer
       (
         "approximate_loops",
         "true approximate loops - false otherwise",
         io::Serialization::GetAgent( &m_ShouldApproximateLoops)
       );

      parameters.AddInitializer
       (
         "use_sans",
         "true use sans implementation - false use saxs implementation",
         io::Serialization::GetAgent( &m_ShouldUseSans)
       );

      parameters.AddOptionalInitializer
      (
        "deuterium_percent",
        "the percent of deuterium in the solvent",
        io::Serialization::GetAgent( &m_DeuteriumExchangeParameter)
      );

      return parameters;
    }
  } // namespace restraint
} // namespace bcl
