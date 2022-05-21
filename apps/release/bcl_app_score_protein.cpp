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
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_vector.h"
#include "assemble/bcl_assemble_aa_sasa_ols.h"
#include "assemble/bcl_assemble_collector_common_aa.h"
#include "assemble/bcl_assemble_quality.h"
#include "assemble/bcl_assemble_quality_batch.h"
#include "assemble/bcl_assemble_sse_geometry_packers.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "biol/bcl_biol_exposure_prediction.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static_and_dynamic.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_extension.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_or.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "contact/bcl_contact_recovery.h"
#include "contact/bcl_contact_statistics.h"
#include "density/bcl_density_simulators.h"
#include "fold/bcl_fold_default_flags.h"
#include "fold/bcl_fold_protocol_em.h"
#include "fold/bcl_fold_protocol_restraint.h"
#include "fold/bcl_fold_protocols.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_binary_function_bind_second.h"
#include "math/bcl_math_roc_curve.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_running_min_max.h"
#include "math/bcl_math_statistics.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "restraint/bcl_restraint_body.h"
#include "restraint/bcl_restraint_contains_body_origin.h"
#include "restraint/bcl_restraint_handler_body.h"
#include "score/bcl_score_aa_neighborhood_exposure.h"
#include "score/bcl_score_aa_neighborhood_exposure_prediction.h"
#include "score/bcl_score_aa_pair_distance_smooth.h"
#include "score/bcl_score_aa_sequence_pair.h"
#include "score/bcl_score_body_assignment.h"
#include "score/bcl_score_density_profile_sse_agreement.h"
#include "score/bcl_score_protein_model_aa_neighborhood.h"
#include "score/bcl_score_protein_model_membrane_topology.h"
#include "score/bcl_score_protein_model_sse_pairs.h"
#include "score/bcl_score_restraint_body_protein_model.h"
#include "score/bcl_score_sse_pair_angle_distance.h"
#include "score/bcl_score_sse_pair_packing.h"
#include "score/bcl_score_strand_pairing.h"
#include "score/bcl_score_symmetry.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "storage/bcl_storage_table.hpp"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinScore
    //! @brief Class for scoring pdbs
    //! @details Class for scoring pdbs according to bcl scoring functions
    //!
    //! @see @link example_app_score_protein.cpp @endlink
    //! @author woetzen, karakam, weinerbe
    //! @date 03/24/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinScore :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! input pdb
      util::ShPtr< command::FlagInterface> m_PDBFileFlag;

      //! input pdb filelist
      util::ShPtr< command::FlagStatic>    m_PDBFileListFlag;
      util::ShPtr< command::Parameter>     m_ListRangeParam;

      //! flag indicating that the list file contains pdbtm xml file containing the membrane and transformation matrix
      //! for each pdb
      util::ShPtr< command::FlagInterface> m_PDBTM_XML_flag;

      //! detailed scoring function output flag
      util::ShPtr< command::FlagInterface> m_DetailedFlag;

      //! flag for inputting a score table instead of calculating scores
      util::ShPtr< command::FlagInterface> m_ScoreTableReadFlag;

      //! flag for outputting scores to table
      util::ShPtr< command::FlagInterface> m_ScoreTableWriteFlag;

      //! flag for calculating enrichments
      util::ShPtr< command::FlagStatic> m_EnrichmentFlag;
      util::ShPtr< command::ParameterInterface> m_EnrichmentGoodModelsFractionParam;
      util::ShPtr< command::ParameterInterface> m_EnrichmentGoodModelsCriteriaCutoffParam;
      util::ShPtr< command::ParameterInterface> m_EnrichmentNumberCrossValidationParam;
      util::ShPtr< command::ParameterInterface> m_EnrichmentCriteriaNameParam;
      util::ShPtr< command::ParameterInterface> m_EnrichmentSortOrderParam;
      util::ShPtr< command::ParameterInterface> m_EnrichmentOutputTableFilenameParam;
      util::ShPtr< command::ParameterInterface> m_EnrichmentPlotFilePrefixParam;

      //! flag to pass file containing math::Vector with weight set to combine score
      util::ShPtr< command::FlagInterface> m_WeightSetFlag;

      //! flag to idealize pdb
      util::ShPtr< command::FlagInterface> m_IdealizeFlag;

      //! flag for contact statistics
      util::ShPtr< command::FlagStatic> m_ContactFlag;

      //! flag for contact statistics
      util::ShPtr< command::FlagStatic> m_DSSPNativeOnlyFlag;

      //! flag for density profile agreement
      util::ShPtr< command::FlagStaticAndDynamic> m_DensityProfileAgreementFlag;
      util::ShPtr< command::ParameterInterface> m_MRCMapParam;
      util::ShPtr< command::ParameterInterface> m_MRCResolutionParam;
      util::ShPtr< command::ParameterInterface> m_SimulatorParameter;

      //! sse prediction methods
      mutable storage::Set< sspred::Method> m_Methods;

      mutable util::ShPtr< biol::Membrane> m_Membrane;

      //! individual scores stored in a table
      mutable storage::Table< double> m_ScoringTable;

      //! weight set
      mutable storage::Table< double> m_Weights;

      //! template/native structure ( for rmsd calculations)
      mutable util::ShPtr< assemble::ProteinModel> m_TemplateModel;

      //! aa pair distance score smooth
      mutable fold::Score e_ScoreAAPairDistanceSmooth;

      mutable score::ProteinModelScoreSum m_ScoreFunction;

      //! density map for density agreement scoring
      mutable util::ShPtr< density::Map> m_DensityMap;

      //! set of atoms that will be used for quality calculations
      mutable storage::Set< biol::AtomType> m_AtomTypes;

      //! qualities that will be calculated
      mutable storage::Set< quality::Measure> m_Qualities;

      //! protein model data to be applied to each protein (i.e. restraints, membrane)
      mutable util::ShPtr< assemble::ProteinModelData> m_ProteinModelData;

      //! the range that is used from the list
      mutable math::Range< size_t> m_ListRange;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinScore();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      ProteinScore *Clone() const
      {
        return new ProteinScore( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief return the original name of the application, as it was used in license files
      //! @return string for the bcl::commons name of that application
      //! This is necessary so that, if release application names change, licenses will continue to work
      std::string GetLicensedName() const
      {
        return "ScoreProtein";
      }

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const
      {
        return "Scores proteins/PDBs according to BCL energy functions";
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      storage::Vector< storage::Triplet< std::string, util::ShPtr< biol::Membrane>, math::TransformationMatrix3D> >
      PDBFilesWithMembraneFromCommandLine() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief returns web text information
      //! @return text (html allowed but not required) that will be displayed on the website
      //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
      const std::string &GetWebText() const;

      //! @brief initialize scores and weights
      void InitializeScores() const;

      //! @brief perform scoring of a single pdb, given by filename - put result in row
      //! @param FILENAME file name of pdb to be scored
      //! @param MEMBRANE the membrane that will be used
      //! @param TRANSFORMATION will be applied to protein before it is scored - relevant for membrane protein models
      //! @param ROW the row results will be written to
      void
      ScoreSinglePDB
      (
        const std::string &FILENAME,
        const biol::Membrane &MEMBRANE,
        const math::TransformationMatrix3D &TRANSFORMATION,
        storage::Row< double> &ROW
      ) const;

      //! @brief check protein validity
      //! @details sometimes protein with unreasonable coordinates are supplied, e.g. output from rosetta loop building
      //! where loops could not have been closed
      //! @param PROTEIN the protein to be scored
      //! @return true, if protein is considered valid
      bool CheckProteinValidity( const assemble::ProteinModel &PROTEIN) const;

      //! @brief determine rank
      //! @return Table with rank of pdb given in command line fore each score in m_Score_table
      storage::Table< size_t> DetermineRank() const;

      //! @brief calculate enrichments
      //! @return Table with same header, but enrichment values in each of the score columns
      storage::Table< double> CalculateEnrichments() const;

      //! @brief calcualtes statistics min, mas, mean and sd for all scores and other stats
      //! @return table with cols as the score and rows for sd, mean, min and max
      storage::Table< double> CalculateStatistics() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }

    private:

      static const ApplicationType ProteinScore_Instance;

    }; // class ScoreProtein

    const ApplicationType ProteinScore::ProteinScore_Instance
    (
      GetAppGroups().AddAppToGroup( new ProteinScore(), GetAppGroups().e_Protein)
    );

    storage::Vector< storage::Triplet< std::string, util::ShPtr< biol::Membrane>, math::TransformationMatrix3D> >
    ProteinScore::PDBFilesWithMembraneFromCommandLine() const
    {
      storage::Vector< storage::Triplet< std::string, util::ShPtr< biol::Membrane>, math::TransformationMatrix3D> > pdb_name_membrane_transformation;

      if( m_ListRangeParam->GetWasSetInCommandLine())
      {
        std::stringstream range_stream( m_ListRangeParam->GetValue());
        std::stringstream error_stream;
        BCL_Assert( m_ListRange.FromStream( range_stream, error_stream), "Error while reading range: " + m_ListRangeParam->GetValue() + " " + error_stream.str());
      }

      //pdbname given as argument
      if( m_PDBFileFlag->GetFlag())
      {
        pdb_name_membrane_transformation.PushBack
        (
          storage::Triplet< std::string, util::ShPtr< biol::Membrane>, math::TransformationMatrix3D>
          (
            m_PDBFileFlag->GetFirstParameter()->GetValue(),
            util::ShPtr< biol::Membrane>
            (
              new biol::Membrane( biol::Membrane::GetCommandLineMembrane())
            ),
            math::TransformationMatrix3D()
          )
        );
      }

      // list of pdbfiles given as argument
      if( m_PDBFileListFlag->GetFlag())
      {
        io::IFStream pdblist;
        io::File::MustOpenIFStream( pdblist, m_PDBFileListFlag->GetFirstParameter()->GetValue());
        size_t count( 0);
        while( !pdblist.eof())
        {
          storage::Triplet< std::string, util::ShPtr< biol::Membrane>, math::TransformationMatrix3D> name_membrane_trans;
          pdblist >> name_membrane_trans.First();

          // sometimes read.eof() does not return false, even if end is reached
          // so it is necessary to check if a name was actually read
          if( name_membrane_trans.First().empty())
          {
            break;
          }

          // if scoring membrane proteins
          // if pdbtm xml files are expected in list
          if( m_PDBTM_XML_flag->GetFlag())
          {
            std::string xml_file;
            pdblist >> xml_file;
            io::IFStream read_xml;
            io::File::MustOpenIFStream( read_xml, xml_file);
            const storage::Pair< biol::Membrane, math::TransformationMatrix3D>
              membrane_transformation
              (
                biol::Membrane::MembraneAndTransformationFromPDBTMXML
                (
                  read_xml,
                  biol::Membrane::GetParameterTransitionThickness()->GetNumericalValue< double>(),
                  biol::Membrane::GetParameterGapThickness()->GetNumericalValue< double>()
                )
              );
            name_membrane_trans.Second() = util::ShPtr< biol::Membrane>
            (
              membrane_transformation.First().Clone()
            );
            name_membrane_trans.Third() = membrane_transformation.Second();
          }
          // score membrane proteins, but xml is not given
          else if( biol::Membrane::GetFlagMembrane()->GetFlag())
          {
            name_membrane_trans.Second() = util::ShPtr< biol::Membrane>
            (
              new biol::Membrane( biol::Membrane::GetCommandLineMembrane())
            );
          }
          // soluble proteins, membrane thickness is 0
          else
          {
            name_membrane_trans.Second() = util::ShPtr< biol::Membrane>( new biol::Membrane( 0.0, 0.0, 0.0));
          }

          if( !m_ListRangeParam->GetWasSetInCommandLine() || m_ListRange.IsWithin( count))
          {
            // store current name and membrane
            pdb_name_membrane_transformation.PushBack( name_membrane_trans);
          }
          ++count;
        }

        io::File::CloseClearFStream( pdblist);
      }

      // end
      return pdb_name_membrane_transformation;
    }

    int ProteinScore::Main() const
    {
      // sspred methods
      if( sspred::Methods::GetFlagReadSSPredictions()->GetFlag())
      {
        m_Methods = sspred::Methods::GetCommandLineMethods();
      }

      const storage::Vector< storage::Triplet< std::string, util::ShPtr< biol::Membrane>, math::TransformationMatrix3D> >
        pdbfiles_membrane( PDBFilesWithMembraneFromCommandLine());

      // initialize the scores
      InitializeScores();

      // initialize weight set
      fold::ScoreWeightSet score_weight_set;

      // set weights from file
      if( m_WeightSetFlag->GetFirstParameter()->GetWasSetInCommandLine())
      {
        // initialize read
        BCL_MessageStd( "reading score weightset from file");
        io::IFStream read;
        io::File::MustOpenIFStream( read, m_WeightSetFlag->GetFirstParameter()->GetValue());
        m_Weights.ReadFormatted( read);
        // close stream
        io::File::CloseClearFStream( read);
        score_weight_set = fold::ScoreWeightSet( m_Weights);
      }
      // from the names in the scoring functions map
      else
      {
        for
        (
          util::Enumerate< util::ShPtr< score::ProteinModel>, fold::Scores>::const_iterator
            itr( fold::GetScores().Begin()), itr_end( fold::GetScores().End());
          itr != itr_end; ++itr
        )
        {
          if( !( **itr)->IsQualityScore())
          {
            score_weight_set.SetWeight( *itr, 1.0);
          }
        }

        // create weights tables with weight one for each of the scores
        m_Weights = score_weight_set.CreateTable();
      }

      // construct the scoring function
      m_ScoreFunction = *( score_weight_set.ConstructScoreSum());

      // if template model flag is given read in the template model
      if( fold::DefaultFlags::GetFlagNativeModel()->GetFlag())
      {
        // initialize read
        io::IFStream read;

        // open stream for native pdb
        io::File::MustOpenIFStream( read, fold::DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetValue());
        pdb::Handler pdb_orig( read);
        io::File::CloseClearFStream( read);

        // initialize sse_min_sizes
        storage::Map< biol::SSType, size_t> sse_min_sizes;
        sse_min_sizes[ biol::GetSSTypes().HELIX] = 0;
        sse_min_sizes[ biol::GetSSTypes().STRAND] = 0;
        sse_min_sizes[ biol::GetSSTypes().COIL] = 0;

        if( m_DSSPNativeOnlyFlag->GetFlag())
        {
          pdb::Factory::GetFlagDSSP()->SetFlag();
        }
        // create model from pdb and store it in the member variable
        m_TemplateModel =
          util::ShPtr< assemble::ProteinModel>
          (
            new assemble::ProteinModel( pdb::Factory().ProteinModelFromPDB( pdb_orig, sse_min_sizes))
          );
        if( m_DSSPNativeOnlyFlag->GetFlag())
        {
          pdb::Factory::GetFlagDSSP()->UnsetFlag();
        }
        // insert the native model into data
        m_ProteinModelData->Insert( assemble::ProteinModelData::e_NativeModel, m_TemplateModel);

        // make a copy of the template model
        util::ShPtr< assemble::ProteinModel> sp_native_filtered_model( m_TemplateModel->HardCopy());
        sp_native_filtered_model->FilterByMinSSESizes( pdb::Factory::GetCommandlineSSETypeMinSizes());
        m_ProteinModelData->Insert( assemble::ProteinModelData::e_NativeFilteredModel, sp_native_filtered_model);

        // read ss predictions
        if( !m_Methods.IsEmpty())
        {
          // try to read all the specified methods in, if it fails warn user
          if
          (
            !sspred::MethodHandler::ReadPredictionsForProteinModel
            (
              m_Methods,
              *m_TemplateModel,
              fold::DefaultFlags::GetFlagReadSequenceDataPath()->GetParameterList()( 1)->GetValue(),
              fold::DefaultFlags::GetFlagReadSequenceDataPath()->GetFirstParameter()->GetValue()
            )
            ||
            !sspred::MethodHandler::ReadPredictionsForProteinModel
            (
              m_Methods,
              *sp_native_filtered_model,
              fold::DefaultFlags::GetFlagReadSequenceDataPath()->GetParameterList()( 1)->GetValue(),
              fold::DefaultFlags::GetFlagReadSequenceDataPath()->GetFirstParameter()->GetValue()
            )
          )
          {
            BCL_MessageCrt( "Can't read all SSMethods for pdb " + fold::DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetValue());
          }
        }

        // read environment predictions
        if( score::AANeighborhoodExposurePrediction::GetFlagScoreExposure()->GetFlag())
        {
          biol::ExposurePrediction::ReadPredictions
          (
            *m_TemplateModel,
            fold::DefaultFlags::GetFlagReadSequenceDataPath()->GetParameterList()( 1)->GetValue(),
            fold::DefaultFlags::GetFlagReadSequenceDataPath()->GetFirstParameter()->GetValue()
          );
          biol::ExposurePrediction::ReadPredictions
          (
            *sp_native_filtered_model,
            fold::DefaultFlags::GetFlagReadSequenceDataPath()->GetParameterList()( 1)->GetValue(),
            fold::DefaultFlags::GetFlagReadSequenceDataPath()->GetFirstParameter()->GetValue()
          );
        }

        // collect the atom types
        m_AtomTypes = biol::AtomTypes::GetCommandLineAtoms();

        // initialize the qualities
        m_Qualities = quality::Measures::GetCommandLineQualityMeasures();

        // if the flag for using native SSE definitions were given
        if
        (
          fold::DefaultFlags::GetFlagUseNativeSSEsAsPool()->GetFlag()
          || !assemble::SSEPool::GetFlagPoolRead()->GetFlag()
        )
        {
          // initialize empty pool
          util::ShPtr< assemble::SSEPool> sp_pool( new assemble::SSEPool());

          // update the pool to have the corresponding SSE definitions from the native model
          sp_pool = util::ShPtr< assemble::SSEPool>
          (
            new assemble::SSEPool
            (
              m_TemplateModel->GetSSEs(),
              true,
              fold::DefaultFlags::GetFlagUseNativeSSEsAsPool()->GetFirstParameter()->GetValue() == "ideal"
            )
          );
          sp_pool->Prune( pdb::Factory::GetCommandlineSSETypeMinSizes());
          m_ProteinModelData->Insert( assemble::ProteinModelData::e_Pool, sp_pool);
        }
      }

      if
      (
        assemble::SSEPool::GetFlagPoolRead()->GetFlag()
        && !fold::DefaultFlags::GetFlagUseNativeSSEsAsPool()->GetFlag()
        && m_TemplateModel.IsDefined()
      )
      {
        // initialize empty pool
        util::ShPtr< assemble::SSEPool> sp_pool( new assemble::SSEPool());

        // open pool file
        io::IFStream read;
        const std::string pool_file( assemble::SSEPool::GetFlagPoolRead()->GetFirstParameter()->GetValue());
        BCL_MessageCrt( "Reading pool from file " + pool_file);
        io::File::MustOpenIFStream( read, pool_file);

        // initialize map to hold the min pool lengths
        storage::Map< biol::SSType, size_t> min_pool_sse_lengths( assemble::SSEPool::GetCommandLineMinSSELengths());

        // read pool
        sp_pool->ReadSSEPool
        (
          read,
          *m_TemplateModel,
          min_pool_sse_lengths[ biol::GetSSTypes().HELIX],
          min_pool_sse_lengths[ biol::GetSSTypes().STRAND]
        );
        io::File::CloseClearFStream( read);

        m_ProteinModelData->Insert( assemble::ProteinModelData::e_Pool, sp_pool);
      }
      // check if scoring table will be read
      if( m_ScoreTableReadFlag->GetFlag())
      {
        BCL_MessageStd
        (
          "read table from file instead of scoring: " + m_ScoreTableReadFlag->GetFirstParameter()->GetValue()
        );
        io::IFStream read;
        io::File::MustOpenIFStream( read, m_ScoreTableReadFlag->GetFirstParameter()->GetValue());
        m_ScoringTable.ReadFormatted( read);
        io::File::CloseClearFStream( read);
      }
      // actual need to score the proteins
      else
      {
        BCL_MessageStd( "scoring all pdbs: " + util::Format()( pdbfiles_membrane.GetSize()));

        // initialize the column names with the default ones
        storage::Vector< std::string> table_header_col_names
        (
          storage::Vector< std::string>::Create
          (
            "nr_aa", "nr_aa_def", "nr_aa_sse", "nr_aa_helix", "nr_aa_strand", "nr_sse", "nr_helix", "nr_strand"
          )
        );

        // if contact stats were requested
        if( m_ContactFlag->GetFlag())
        {
          table_header_col_names.PushBack( "per_co_short");
          table_header_col_names.PushBack( "per_co_mid");
          table_header_col_names.PushBack( "per_co_long");
        }

        // add quality names
        table_header_col_names.Append( assemble::QualityBatch::ColumnNamesFromQualities( m_Qualities));

        // add score names and the sum
        table_header_col_names.Append( m_Weights.GetHeader());
        table_header_col_names.PushBack( "sum");

        // initialize the scoring table with the table header constructed from the column names
        m_ScoringTable = storage::Table< double>( storage::TableHeader( table_header_col_names));

        // iterate over the pdbs
        for
        (
          storage::Vector< storage::Triplet< std::string, util::ShPtr< biol::Membrane>, math::TransformationMatrix3D> >::const_iterator
            filename_itr( pdbfiles_membrane.Begin()), filename_itr_end( pdbfiles_membrane.End());
          filename_itr != filename_itr_end;
          ++filename_itr
        )
        {
          // get the row where the results should be written to
          storage::Row< double> &current_result( m_ScoringTable.InsertRow( filename_itr->First(), true));

          // score the PDB, while also writing the results to the corresponding row
          ScoreSinglePDB( filename_itr->First(), *filename_itr->Second(), filename_itr->Third(), current_result);
        }

        // only write if the flag was set
        if( m_ScoreTableWriteFlag->GetFlag())
        {
          io::OFStream write;
          io::File::MustOpenOFStream( write, m_ScoreTableWriteFlag->GetFirstParameter()->GetValue());
          m_ScoringTable.WriteFormatted( write, util::Format().W( 10));
          io::File::CloseClearFStream( write);
        }
        // print to screen
        else
        {
          m_ScoringTable.WriteFormatted( util::GetLogger(), util::Format().W( 10));
        }
      }

      // if enrichment should be calculated
      if( m_EnrichmentFlag->GetFlag())
      {
        const storage::Table< double> enrichments( CalculateEnrichments());
        BCL_MessageStd( "calculated enrichments");
        const std::string output_table_filename( m_EnrichmentOutputTableFilenameParam->GetValue());
        io::OFStream write;
        io::File::MustOpenOFStream( write, output_table_filename);
        enrichments.WriteFormatted( write, util::Format().W( 12));
        io::File::CloseClearFStream( write);
        enrichments.GetTransposedTable().WriteFormatted( util::GetLogger(), util::Format().W( 12));
      }

      // if m_Rank_flag was given
      if( m_TemplateModel.IsDefined())
      {
        const storage::Table< size_t> rank_table( DetermineRank());
        BCL_MessageStd( "calculated ranks for pdb: " + fold::DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetValue());
        rank_table.WriteFormatted( util::GetLogger(), util::Format().W( 12));
      }

      // get the statistics table
      const storage::Table< double> stats_table( CalculateStatistics());

      // write statistics
      BCL_MessageStd( "statistics table");
      stats_table.WriteFormatted( util::GetLogger(), util::Format().W( 12));

      // write transposed statistics
      BCL_MessageStd( "transposed statistics table");
      stats_table.GetTransposedTable().WriteFormatted( util::GetLogger(), util::Format().W( 12));

      return 0;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProteinScore::GetReadMe() const
    {
      // create the readme text
      static const std::string s_readme_text
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This application scores proteins using the BCL protein scoring function, which is based on knowledge based"
        "potential functions.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::Score?\n"
        "BCL::Score is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part of "
        "a larger library of applications called BCL::Commons.\n"
        "BCL::Score rates given protein structures according to the likelihood of this conformation appearing in nature."
        "Therefore BCL::Score uses knowledge based potential functions based on statistics collected from known protein"
        "structures. These potential include amino acid environments, amino acid pair distance, SSE pairing/packing, "
        "loop lengths as well as membrane specific properties.\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::Score.\n"
        "When using BCL::Score in a publication, please cite the publications describing the application's "
        "development"
        "\n"
        "N. Woetzel, M. Karakaş, R. Staritzbichler, R. Müller, B. E. Weiner, and J. Meiler, BCL::Score--knowledge "
        "based energy potentials for ranking protein models represented by idealized secondary structure elements., "
        "PLoS One, vol. 7, no. 11, p. e49242, Jan. 2012.\n"
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL::Score.\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator()
      );

      return s_readme_text;
    }

    //! @brief returns web text information
    //! @return text (html allowed but not required) that will be displayed on the website
    //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
    const std::string &ProteinScore::GetWebText() const
    {
      // create a static string to hold readme information
      static const std::string s_web_text
      (
        "BCL::Score: Scores protein structures according to the likelihood of a given confirmation appearing in nature\n\n"
        "This application uses a knowledge-based energy function built from statistics collected from known protein "
        "structures.  These statistics are used to build potentials including amino acid environments, amino acid pair "
        "distance, SSE pairing/packing, loop lengths, radius of gyration, contact order, and SSE prediction agreement.\n\n"
        "The aim of BCL::Score is to identify native-like protein topology for arrangements of secondary structure "
        "elements only. The rationale for this approach is the hypothesis that interactions between SSEs define the core"
        "of the protein structure and are the major contributor to the stability of the protein fold.  Further, these "
        "stabilizing interactions can be accurately predicted as flexibility is reduced in the backbone of SSEs without "
        "side chains and loop regions.  Because of this loss in noise, we expect a smoothened energy landscape that is  "
        "readily searched. \n"
        "!score1.png!\n"
        "Fig. 1: Amino acid neighbor count environment potential.  A shows the transition function that is used between "
        "the lower and upper threshold in which the weight for the neighbor being considered drops from 1 ( 4 angstrom) "
        "to 0 ( 11.4 angstrom) using half of a cosine function.  B shows the neighbor count energy potential for all 20 "
        "amino acids with their three letter code. \n\n"
        "!score2.png!\n"
        "Fig. 2: Amino acid pair distance potentials. In A the idealized structure of 1ubi with C-beta and H-alpha 2 "
        "atoms is shown with the distances between ILE 32 and LEU 56 (4.7 angstrom) and between LYS 11 and GLU 34 "
        "(8.3 angstrom). B shows selected amino acid pair distance potentials for Trp-Trp as an example for pie stacking "
        "interaction, ILE-LEU as an example for vdW apolar interaction, ARG-GLU as an example for Coulomb attraction "
        " and ARG-LYS as an example for Coulomb repulsion. \n\n"
        "!score3.png!\n"
        "Fig. 3: Loop closure potential. A describes two beta-strands connected by a loop characterized by the Euclidean "
        "distance between the tow ends and the number of residues in the loop connecting these two ends.  B describes the "
        "derived energy potential, where the energy is a function of the number of residues in the loop and the Euclidean "
        "distance between the ends of the main axes. \n\n"
        "!score4.png!\n"
        "Fig. 4: SSE Fragment packing.  SSE fragments are shown with their geometric packing descriptors.  A alpha1 and "
        "alpha2 are orthogonal, if the shortest connection between the main axes is orthogonal.  B connection is not "
        "orthogonal, since the minimal interface length m cannot be achieved.  C theta is the twist angle around the "
        "shortest connection - which is equivalent to the dihedral angle between main axis 1 - shortest connection - "
        "main axis 2.  D omega is hte offset form the optimal expected position for a helix-strand interaction, if it "
        "is 0 degrees, the helix is on top of the strand, if it is 90 degrees, the helix would interact with the backbone "
        "of the strand.  Omega1 and omega2 are the offsets for a strand-strand packing - for omegas close to 90 degrees, "
        "it is a strand backbone pairing interaction dominated  by hydrogen bond interaction within a sheet, if they are "
        "close to 0 degrees, it is dominated by side chain interactions like seen in sheet-sandwiches.  E every SSE is "
        "represented as multiple fragments and the SSE interaction is described by the list of all fragment interactions, "
        "leaving out additional fragments of the longer SSE with sub-optimal packing (bottom grey helix fragment). \n\n"
        "!score5.png!\n"
        "Fig. 5: Strand pairing and SSE packing potential.  Shown are all secondary structure element packing potentials "
        "with their schematic shortest connections, twist angle and their derived potentials.  A shows the beta-strand-beta-strand "
        "pairing potential with prominent distance of 4.75 angstrom and angles of -15 degrees and 165 degrees.  B shows "
        "the alpha-Helix-alpha-Helix packing with preferred packing distance of 10 angstrom and the preferred parallel "
        "angle of -45 degrees and the anti-parallel packing of 135 degrees.  C shows the Beta-Sheet-Beta-Sheet packing "
        "potential with a preferred distance 10 angstroms and angles of -30 degrees and 150 degrees. D shows the alpha-"
        "Helix-Beta-Sheet packing with its packing distance around 10 angstrom and an anti-parallel angle of 150-180 "
        "degrees. \n\n"
        "!score6.png!\n"
        "Fig. 6: Contact order and square radius of gyration potential.  A Fold complexity is represented by the contact "
        "order potential.  The potential is given as the likelihood to observe a contact order to number of residues ratio "
        "in the model.  B Statistics for the square radius of gyration over the number of residues were directly collected "
        "in a histogram and converted into a potential. \n\n"
      );
      return s_web_text;
    }

    //! @brief initialize scores and weights
    void ProteinScore::InitializeScores() const
    {
      // initialize only once
      e_ScoreAAPairDistanceSmooth = fold::Score( score::AAPairDistanceSmooth::GetDefaultScheme());
      if( e_ScoreAAPairDistanceSmooth.IsDefined())
      {
        return;
      }

      // AAPairDistanceSmooth
      e_ScoreAAPairDistanceSmooth = fold::GetScores().AddScore
      (
        util::ShPtr< score::ProteinModel>
        (
          new score::ProteinModelSSEPairs
          (
            util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >
            (
              new score::AASequencePair( score::AAPairDistanceSmooth(), false)
            ),
            false
          )
        )
      );

      // initialize all scores from all protocols
      for( fold::Protocols::iterator itr( fold::GetProtocols().Begin()), itr_end( fold::GetProtocols().End()); itr != itr_end; ++itr)
      {
        ( **itr)->InitializeScores();
      }

      // exposures
      fold::GetScores().AddScore
      (
        util::ShPtr< score::ProteinModel>
        (
          new score::ProteinModelAANeighborhood
          (
            util::CloneToShPtr( score::AANeighborhoodExposure( assemble::AANeighborVector()))
          )
        )
      );
      fold::GetScores().AddScore
      (
        util::ShPtr< score::ProteinModel>
        (
          new score::ProteinModelAANeighborhood
          (
            util::CloneToShPtr( score::AANeighborhoodExposure( assemble::AASasaOLS()))
          )
        )
      );

      // secondary structure elements
      fold::GetScores().AddScore
      (
        util::ShPtr< score::ProteinModel>
        (
          new score::ProteinModelSSEPairs
          (
            util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >
            (
              new score::SSEPairAngleDistance( score::SSEPairPacking())
            ),
            false
          )
        )
      );
      fold::GetScores().AddScore
      (
        util::ShPtr< score::ProteinModel>
        (
          new score::ProteinModelSSEPairs
          (
            util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >
            (
              new score::SSEPairAngleDistance( score::StrandPairing())
            ),
            false
          )
        )
      );

      // density profile
      if( m_DensityProfileAgreementFlag->GetFlag())
      {
        // check that resolution and map were given
        const double mrc_resolution( m_MRCResolutionParam->GetNumericalValue< double>());

        // create HandlerAtomDistanceAssigned "handler" as method for determining if a restraint body is occupied
        restraint::HandlerBody handler
        (
          util::ShPtr< util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool> >
          (
            new restraint::ContainsBodyOrigin()
          )
        );

        // create string "restraint_filename" from the native structure
        const std::string restraint_filename( fold::DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetValue());

        // create stream to restraint file
        io::IFStream read;
        io::File::MustOpenIFStream( read, restraint_filename);
        BCL_MessageStd( "read restraint file: " + restraint_filename);

        // create restraints
        const util::ShPtr< restraint::Body> restraints( handler.CreateRestraintsBody( read).FirstElement());
        io::File::CloseClearFStream( read);

        const util::ShPtr< density::SimulateInterface> simulator
        (
          density::GetSimulators().CreateSimulator
          (
            density::Simulator( m_SimulatorParameter->GetValue()), m_DensityMap->GetCellWidth(), mrc_resolution
          )
        );

        // iterate over dynamic parameters
        const util::SiPtrVector< const command::ParameterInterface> dynamic_parameters( m_DensityProfileAgreementFlag->GetDynamicParameterList());

        // iterate over profile types
        for
        (
          util::SiPtrVector< const command::ParameterInterface>::const_iterator
            itr( dynamic_parameters.Begin()), itr_end( dynamic_parameters.End());
          itr != itr_end;
          ++itr
        )
        {
          // construct from simulator, density and restraint
          util::ShPtr< score::DensityProfileSSEAgreement> score_profile_agreement
          (
            new score::DensityProfileSSEAgreement
            (
              simulator,
              *m_DensityMap,
              *restraints,
              score::DensityProfileSSEAgreement::Profile1DTypeEnum( ( *itr)->GetValue())
            )
          );

          // create protein model score
          const util::ShPtr< score::RestraintBodyProteinModel> score_model
          (
            new score::RestraintBodyProteinModel
            (
              util::CloneToShPtr
              (
                util::ShPtrVector< restraint::Body>( 1, restraints)
              ),
              // create an assignment score for angle_1D
              score::BodyAssignment( score_profile_agreement)
            )
          );

          fold::GetScores().AddScore( score_model);
        }
      }

      // if membrane flag is passed
      if( biol::Membrane::GetFlagMembrane()->GetFlag())
      {
        // update the protein model data
        m_ProteinModelData->Insert( assemble::ProteinModelData::e_Membrane, m_Membrane);
      }
    }

    //! @brief perform scoring of a single pdb, given by filename
    void
    ProteinScore::ScoreSinglePDB
    (
      const std::string &FILENAME,
      const biol::Membrane &MEMBRANE,
      const math::TransformationMatrix3D &TRANSFORMATION,
      storage::Row< double> &ROW
    ) const
    {
      // update the membrane
      *m_Membrane = MEMBRANE;

      //instantiate pdb
      io::IFStream read;
      io::File::MustOpenIFStream( read, FILENAME);
      BCL_MessageTop( "scoring pdb: " + FILENAME);
      pdb::Handler pdb( read);

      io::File::CloseClearFStream( read);

      // create model from pdb
      assemble::ProteinModel model( pdb::Factory().ProteinModelFromPDB( pdb));

      // check for read in membrane, use that instead if it was found
      util::ShPtr< assemble::ProteinModelData> sp_data( model.GetProteinModelData());
      const util::ShPtr< biol::Membrane> sp_membrane( sp_data->GetData( assemble::ProteinModelData::e_Membrane));
      if( sp_membrane.IsDefined() && !m_PDBTM_XML_flag->GetFlag())
      {
        m_Membrane = sp_membrane;
      }

      // update protein model data
      sp_data->Insert( *m_ProteinModelData, true);

      // if there is protein model data
      if( !sp_data->IsEmpty())
      {
        // add it to the model
        model.SetProteinModelData( sp_data);
      }

      // if the model has no chains then skip
      if( model.GetNumberOfChains() == 0)
      {
        BCL_MessageCrt( "pdb does not contain any chain");
        return;
      }

      // if protein is not a valid structure, has long disconnects between CA residue then skip
      if( !CheckProteinValidity( model))
      {
        BCL_MessageCrt( "protein " + FILENAME + " is not considered a valid protein! Skipping.");
        return;
      }

      // transform the protein model according to given transformation if membrane protein scoring is asked for
      if( biol::Membrane::GetFlagMembrane()->GetFlag() || m_PDBTM_XML_flag->GetFlag())
      {
        model.Transform( TRANSFORMATION);
      }

      // only one sse
      if
      (
        !pdb::Factory::GetFlagDSSP()->GetFlag() // if DSSP found no SSEs, lets assume that it's right
        &&
        (
             model.GetNumberSSEs() == 0
          || ( model.GetNumberSSEs() == 1 && model.GetSSEs().FirstElement()->GetType() == biol::GetSSTypes().COIL)
        )
      )
      {
        // new protein model
        assemble::ProteinModel new_model;
        BCL_MessageStd( "using phi psi backbone conformation to assign SSEs in model");
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator
            chain_itr( model.GetChains().Begin()), chain_itr_end( model.GetChains().End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          util::ShPtr< assemble::Chain> new_chain
          (
            new assemble::Chain( assemble::ConstructChainWithSSEsFromConformation( ( *chain_itr)->GetSequence()))
          );
          new_model.Insert( new_chain);
          new_chain->Join( biol::GetSSTypes().COIL, true);
        }
        BCL_MessageStd
        (
          "found " + util::Format()( new_model.GetNumberSSEs()) +
          " sses using ramachandran phi psi angles for identification"
        );

        // update the model
        new_model.SetProteinModelData( model.GetProteinModelData());
        model = new_model;

        // filter the small SSEs
        //model.AddLoops( true, true);
        if( pdb::Factory::GetFlagMinSSESize()->GetFlag())
        {
          model.FilterByMinSSESizes( pdb::Factory::GetCommandlineSSETypeMinSizes());
        }
        BCL_MessageStd( "number of sses in model: " + util::Format()( model.GetNumberSSEs()));
      }

      // if template model is defined
      if( m_TemplateModel.IsDefined())
      {
        // if the model should be superimposed prior to scoring (would change location in membrane, for example)
        if( quality::SuperimposeMeasures::GetFlagSuperimposeMeasure()->GetFlag())
        {
          // construct the superimpose measure
          const quality::SuperimposeMeasure superimpose
          (
            quality::SuperimposeMeasures::GetFlagSuperimposeMeasure()->GetFirstParameter()->GetValue()
          );

          // perform the superimposition
          if( superimpose.IsDefined() && superimpose != quality::GetSuperimposeMeasures().e_NoSuperimpose)
          {
            assemble::Quality::SuperimposeModel( superimpose, model, m_AtomTypes);
          }
        }
      }
      //pdbcode - removing .pdb
      storage::VectorND< 2, std::string> path_code( io::File::SplitToPathAndFileName( FILENAME));
      path_code.Second() = io::File::RemoveFullExtension( path_code.Second());

      // read ss predictions
      if( !m_Methods.IsEmpty())
      {
        // try to read all the specified methods in, if it fails warn user
        if
        (
          !sspred::MethodHandler::ReadPredictionsForProteinModel( m_Methods, model, path_code.Second(), path_code.First())
          &&
          (
            !fold::DefaultFlags::GetFlagReadSequenceDataPath()->GetFlag()
            ||
            !sspred::MethodHandler::ReadPredictionsForProteinModel
            (
              m_Methods,
              model,
              fold::DefaultFlags::GetFlagReadSequenceDataPath()->GetParameterList()( 1)->GetValue(),
              fold::DefaultFlags::GetFlagReadSequenceDataPath()->GetFirstParameter()->GetValue()
            )
          )
        )
        {
          BCL_MessageCrt( "Can't read all SSMethods for pdb " + FILENAME);
        }
      }
      // read environment predictions
      if( score::AANeighborhoodExposurePrediction::GetFlagScoreExposure()->GetFlag())
      {
        biol::ExposurePrediction::ReadPredictions
        (
          model,
          fold::DefaultFlags::GetFlagReadSequenceDataPath()->GetParameterList()( 1)->GetValue(),
          fold::DefaultFlags::GetFlagReadSequenceDataPath()->GetFirstParameter()->GetValue()
        );
      }

    ///////////
    // stats //
    ///////////

      // initialize variables
      size_t nr_residues( 0);
      // get the sequences
      const util::SiPtrVector< const biol::AASequence> sequences( model.GetSequences());
      // iterate over sequences
      for
      (
        util::SiPtrVector< const biol::AASequence>::const_iterator
          seq_itr( sequences.Begin()), seq_itr_end( sequences.End());
        seq_itr != seq_itr_end;
        ++seq_itr
      )
      {
        nr_residues += ( *seq_itr)->GetSize();
      }

      ROW[ "nr_aa"] = nr_residues;
      ROW[ "nr_aa_def"] = assemble::CollectorCommonAA::CollectDefinedAAsInSSEs( model).GetSize();
      size_t nr_sse_residues( 0);
      size_t nr_helix_residues( 0);
      size_t nr_strand_residues( 0);

      // get all the SSEs
      const util::SiPtrVector< const assemble::SSE> sses( model.GetSSEs());
      // iterate over them
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( sses.Begin()), sse_itr_end( sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // store this type
        const biol::SSType this_type( ( *sse_itr)->GetType());
        const size_t this_size( ( *sse_itr)->GetSize());

        // if it's a helix
        if( this_type == biol::GetSSTypes().HELIX)
        {
          // increment helix and SSE amino acid counts
          nr_helix_residues += this_size;
          nr_sse_residues += this_size;
        }
        // if strand
        else if( this_type == biol::GetSSTypes().STRAND)
        {
          // increment strand and SSE amino acid counts
          nr_strand_residues += this_size;
          nr_sse_residues += this_size;
        }
      }
      ROW[ "nr_aa_sse"] = nr_sse_residues;
      ROW[ "nr_aa_helix"] = nr_helix_residues;
      ROW[ "nr_aa_strand"] = nr_strand_residues;

      // statistics for number SSEs
      const size_t nr_helix( model.GetNumberSSE( biol::GetSSTypes().HELIX));
      const size_t nr_strand( model.GetNumberSSE( biol::GetSSTypes().STRAND));
      ROW[ "nr_sse"]    = nr_helix + nr_strand;
      ROW[ "nr_helix"]  = nr_helix;
      ROW[ "nr_strand"] = nr_strand;

      // idealize sses in model
      if( m_IdealizeFlag->GetFlag())
      {
        model.SetToIdealConformation();
      }

      // if contact stats were requested
      if( m_ContactFlag->GetFlag())
      {
        // calculate contact information
        const double ratio_contacts_short( contact::Statistics( contact::Statistics::e_RatioContactsShort, true)( model));
        const double ratio_contacts_mid( contact::Statistics( contact::Statistics::e_RatioContactsMid, true)( model));
        const double ratio_contacts_long( contact::Statistics( contact::Statistics::e_RatioContactsLong, true)( model));

        // update values
        ROW[ "per_co_short"] = ratio_contacts_short;
        ROW[ "per_co_mid"] = ratio_contacts_mid;
        ROW[ "per_co_long"] = ratio_contacts_long;
      }

      // if detailed flag is given
      if( m_DetailedFlag->GetFlag())
      {
        // iterate over all the scores
        m_ScoreFunction.WriteDetailedSchemeAndValues( model, util::GetLogger());
      }

      // actual scores
      const storage::Table< double> scores( m_ScoreFunction.CreateValueTableHorizontal( model));

      // insert all the scores into the ROW
      for
      (
        storage::TableHeader::const_iterator itr( scores.GetHeader().Begin()), itr_end( scores.GetHeader().End());
        itr != itr_end;
        ++itr
      )
      {
        // if their is no such col in the scoring table, just skip
        if( !ROW.GetHeader().HasColumn( *itr))
        {
          continue;
        }

        // copy value from scores to the scoring table
        ROW[ *itr] = scores[ "weighted_value"][ *itr];
      }

      // copy weighted sum
      ROW[ "sum"] = scores[ "weighted_value"][ "sum"];

      // calculate the qualities
      // construct a quality batch and get the list of quality names and the values
      const storage::List< storage::Pair< std::string, double> > quality_values
      (
        assemble::QualityBatch( m_Qualities, m_AtomTypes)( model)
      );

      // iterate over the list of quality names and values
      for
      (
        storage::List< storage::Pair< std::string, double> >::const_iterator
          quality_itr( quality_values.Begin()), quality_itr_end( quality_values.End());
        quality_itr != quality_itr_end; ++quality_itr
      )
      {
        // insert into the ROW
        ROW[ quality_itr->First()] = quality_itr->Second();
      }

      // end
      return;
    }

    //! @brief check protein validity
    //! @details sometimes protein with unreasonable coordinates are supplied, e.g. output from rosetta loop building
    //! where loops could not have been closed
    //! @param PROTEIN the protein to be checked
    //! @return true, if protein is considered valid
    bool ProteinScore::CheckProteinValidity( const assemble::ProteinModel &PROTEIN) const
    {
      // iterate over all chains
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN.GetChains().Begin()), chain_itr_end( PROTEIN.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // get all CA coordinates
        static const storage::Set< biol::AtomType> s_ca_atom_type( biol::GetAtomTypes().CA);
        const util::SiPtrVector< const linal::Vector3D> ca_coordinates( ( *chain_itr)->GetAtomCoordinates( s_ca_atom_type));

        // iterate over coordinates and return false of two defined coordinates are more than 50A apart
        for
        (
          util::SiPtrVector< const linal::Vector3D>::const_iterator
            itr1( ca_coordinates.Begin()), itr_end( ca_coordinates.End());
          itr1 != itr_end;
          ++itr1
        )
        {
          // defined coordinate
          if( !( *itr1)->IsDefined())
          {
            continue;
          }

          const util::SiPtrVector< const linal::Vector3D>::const_iterator itr2( itr1 + 1);
          if( itr2 == itr_end)
          {
            break;
          }
          if( !( *itr2)->IsDefined())
          {
            // skip this pair iterator
            itr1 = itr2;
          }

          // calculate distance
          const double distance( linal::Distance( **itr1, **itr2));

          // unreasonable distance
          static const double s_unreasonable_distance( 100.0);
          if( distance > s_unreasonable_distance)
          {
            return false;
          }
        }
      }

      // end
      return true;
    }

    //! @brief determine rank
    //! @return Table with rank of pdb given in command line fore each score in m_Score_table
    storage::Table< size_t> ProteinScore::DetermineRank
    (
    ) const
    {
      // use same header for the rank table
      // each column will contain the rank of the protein model given with the template model within that score
      storage::Table< size_t> rank_table( m_ScoringTable.GetHeader());

      if( !m_ScoringTable.HasRow( fold::DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetValue()))
      {
        storage::Row< double>
          &current_result( m_ScoringTable.InsertRow( fold::DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetValue()));
        ScoreSinglePDB
        (
          fold::DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetValue(),
          biol::Membrane::GetCommandLineMembrane(),
          math::TransformationMatrix3D(),
          current_result
        );
      }

      // add scores and rmsd to score table
      const storage::Row< double>
        &current_result( m_ScoringTable[ fold::DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetValue()]);
      BCL_MessageStd( "scores of pdb to be ranked:");
      current_result.WriteFormatted( util::GetLogger());

      // insert rank row
      storage::Row< size_t> &current_rank_row( rank_table.InsertRow( fold::DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetValue()));

      // acquire comparison function
      const math::Comparisons< double>::Comparison comp_function( m_EnrichmentSortOrderParam->GetValue());

      // iterate over header and sort column by the score
      for
      (
        storage::TableHeader::const_iterator
          header_itr( rank_table.GetHeader().Begin()), header_itr_end( rank_table.GetHeader().End());
        header_itr != header_itr_end;
        ++header_itr
      )
      {
        m_ScoringTable.SortByColumn( *header_itr, **comp_function);
        current_rank_row[ *header_itr] = m_ScoringTable.RowIndex( fold::DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetValue());
      }

      // end
      return rank_table;
    }

    //! @brief calculate enrichments
    //! @return Table with same header, but enrichment values in each of the score columns
    storage::Table< double> ProteinScore::CalculateEnrichments
    (
    ) const
    {
      // find first iterator that is larger than given rmsd threshold
      const std::string criteria_name( m_EnrichmentCriteriaNameParam->GetValue());
      const double criteria_cutoff( m_EnrichmentGoodModelsCriteriaCutoffParam->GetNumericalValue< double>());
      const size_t criteria_index( m_ScoringTable.GetHeader()[ m_EnrichmentCriteriaNameParam->GetValue()]);

      // sort by given criteria
      const math::Comparisons< double>::Comparison comp_function( m_EnrichmentSortOrderParam->GetValue());
      const math::BinaryFunctionBindSecond< double, double, bool> classify_function( *comp_function, criteria_cutoff);
      m_ScoringTable.SortByColumn( criteria_name, **comp_function);

      // create new table header with adjacent standard deviation
      storage::TableHeader table_header;
      for
      (
        storage::TableHeader::const_iterator
          itr( m_ScoringTable.GetHeader().Begin()), itr_end( m_ScoringTable.GetHeader().End());
        itr != itr_end;
        ++itr
      )
      {
        table_header.PushBack( *itr);
        table_header.PushBack( *itr + "_sd");
      }

      // initialize table with enrichments
      storage::Table< double> enrichments( table_header);

      storage::Table< double>::const_iterator row_itr_split( m_ScoringTable.Begin()), row_itr_end( m_ScoringTable.End());
      while
      (
        row_itr_split != row_itr_end &&
        classify_function( row_itr_split->Second()( criteria_index))
      )
      {
        ++row_itr_split;
      }

      // check that item was found
      if( row_itr_split == row_itr_end)
      {
        BCL_MessageCrt
        (
          "could only find any structures with " + criteria_name + " " + comp_function.GetName() +
          " than given criteria cutoff: " + util::Format()( criteria_cutoff)
        );

        // end
        return enrichments;
      }

      // check that item was found
      if( row_itr_split == m_ScoringTable.Begin())
      {
        BCL_MessageCrt
        (
          "could not find any structures with " + criteria_name + " " + comp_function.GetName() +
          " than given criteria cutoff: " + util::Format()( criteria_cutoff)
        );

        // end
        return enrichments;
      }

      const storage::Table< double> good_models( m_ScoringTable.Begin(), row_itr_split);
      const size_t nr_good_models( good_models.GetSize());
      const storage::Table< double> bad_models( row_itr_split, row_itr_end);
      const size_t nr_bad_models( bad_models.GetSize());
      const double target_fraction( m_EnrichmentGoodModelsFractionParam->GetNumericalValue< double>());

      BCL_MessageStd( "according to criteria, nr good models: " + util::Format()( nr_good_models));
      BCL_MessageStd( "according to criteria, nr bad models: " + util::Format()( nr_bad_models));

      storage::List< storage::Table< double> > cross_validation_list;

      // store number of cross validations
      const size_t nr_cross_validations( m_EnrichmentNumberCrossValidationParam->GetNumericalValue< size_t>());

      // too many bad models
      if( double( nr_good_models) / double( nr_good_models + nr_bad_models) < target_fraction)
      {
        const size_t target_nr_bad_models( size_t( double( nr_good_models) / target_fraction) - nr_good_models);

        cross_validation_list = bad_models.CreateDifferentSubTables( nr_cross_validations, target_nr_bad_models);

        // append good models to each good model in cross validation list
        for
        (
          storage::List< storage::Table< double> >::iterator
            itr( cross_validation_list.Begin()), itr_end( cross_validation_list.End());
          itr != itr_end;
          ++itr
        )
        {
          itr->Append( good_models);
        }
      }
      // too many good models
      else
      {
        // the number of good models we should have after removing the excess ones
        size_t target_nr_good_models( size_t( double( nr_bad_models) / ( 1.0 - target_fraction)) - nr_bad_models);

        // make sure that target_nr_good models is less than nr_good_models - (number_cross_validation -1)
        // otherwise it wouldn't be able to make the specified number of distinct subtables
        if( nr_good_models - target_nr_good_models < nr_cross_validations - 1)
        {
          // update number true entries to have a smaller size
          target_nr_good_models = nr_good_models - nr_cross_validations - 1;
        }

        // initialize the good models lists
        cross_validation_list = good_models.CreateDifferentSubTables( nr_cross_validations, target_nr_good_models);

        // append bad models to each good model in cross validation list
        for
        (
          storage::List< storage::Table< double> >::iterator
            itr( cross_validation_list.Begin()), itr_end( cross_validation_list.End());
          itr != itr_end;
          ++itr
        )
        {
          itr->Append( bad_models);
        }
      }

      storage::Row< double>
        &current_row
         (
           enrichments.InsertRow
           (
             m_EnrichmentGoodModelsCriteriaCutoffParam->GetValue() + "_" + util::Format()( nr_good_models)
           )
         );

      // calculate enrichments for each score
      for
      (
        storage::TableHeader::const_iterator
          itr_score( m_ScoringTable.GetHeader().Begin()),
          itr_score_end( m_ScoringTable.GetHeader().End());
        itr_score != itr_score_end;
        ++itr_score
      )
      {
        linal::Vector< double> each_enrichment( cross_validation_list.GetSize(), double( 0));
        double *ptr( each_enrichment.Begin());
        // calculate enrichment for each item in cross validation set
        for
        (
          storage::List< storage::Table< double> >::iterator
            itr( cross_validation_list.Begin()), itr_end( cross_validation_list.End());
          itr != itr_end;
          ++itr, ++ptr
        )
        {
          itr->SortByColumn( *itr_score);
          const math::ROCCurve
            roc_curve
            (
              itr->ExtractDataPairs( *itr_score, m_EnrichmentCriteriaNameParam->GetValue()),
              classify_function
            );

          *ptr = roc_curve.ContingencyMatrixFraction( target_fraction).GetEnrichment();

          if( m_EnrichmentPlotFilePrefixParam->GetWasSetInCommandLine())
          {
            // filename for plot file
            const std::string plot_file_name
            (
              m_EnrichmentPlotFilePrefixParam->GetValue() + util::Format()( ptr - each_enrichment.Begin()) +
              ( *itr_score) + ".plot"
            );
            io::OFStream write;
            io::File::MustOpenOFStream( write, plot_file_name);
            roc_curve.WriteRatePlottingTable( write);
            io::File::CloseClearFStream( write);
          }
        }

        current_row[ *itr_score] = math::Statistics::Mean( each_enrichment.Begin(), each_enrichment.End());
        current_row[ *itr_score + "_sd"]
                     = math::Statistics::StandardDeviation( each_enrichment.Begin(), each_enrichment.End());
      }

      // end
      return enrichments;
    }

    //! @brief calculates statistics min, max, mean and sd for all scores and other stats
    //! @return table with cols as the score and rows for sd, mean, min and max
    storage::Table< double> ProteinScore::CalculateStatistics() const
    {
      // collects the mins and max
      storage::Map< std::string, math::RunningMinMax< double> > min_max_collector;
      storage::Map< std::string, math::RunningAverageSD< double> > mean_sd_collector;

      // iterate over all rows in the scoring table
      for
      (
        storage::Table< double>::const_iterator row_itr( m_ScoringTable.Begin()), row_itr_end( m_ScoringTable.End());
        row_itr != row_itr_end;
        ++row_itr
      )
      {
        // iterate over all the columns
        for
        (
          storage::TableHeader::const_iterator
            col_itr( m_ScoringTable.GetHeader().Begin()), col_itr_end( m_ScoringTable.GetHeader().End());
          col_itr != col_itr_end;
          ++col_itr
        )
        {
          const double &value( row_itr->Second()[ *col_itr]);
          min_max_collector[ *col_itr] += value;
          mean_sd_collector[ *col_itr] +=  value;
        }
      }

      // create table with the stats
      storage::Table< double> stats( m_ScoringTable.GetHeader());
      stats.InsertRow( "min",  storage::Vector< double>( m_ScoringTable.GetHeader().GetSize(), double( 0)));
      stats.InsertRow( "max",  storage::Vector< double>( m_ScoringTable.GetHeader().GetSize(), double( 0)));
      stats.InsertRow( "mean", storage::Vector< double>( m_ScoringTable.GetHeader().GetSize(), double( 0)));
      stats.InsertRow( "sd",   storage::Vector< double>( m_ScoringTable.GetHeader().GetSize(), double( 0)));

      // iterate over all the columns
      for
      (
        storage::TableHeader::const_iterator
          col_itr( m_ScoringTable.GetHeader().Begin()), col_itr_end( m_ScoringTable.GetHeader().End());
        col_itr != col_itr_end;
        ++col_itr
      )
      {
        // cast to proper pointer type
        const math::RunningMinMax< double> &min_max_collectorb( min_max_collector[ *col_itr]);
        stats[ "min"][ *col_itr] = min_max_collectorb.GetMin();
        stats[ "max"][ *col_itr] = min_max_collectorb.GetMax();
        const math::RunningAverageSD< double> &mean_sd_collectorb( mean_sd_collector[ *col_itr]);
        stats[ "mean"][ *col_itr] = mean_sd_collectorb.GetAverage();
        stats[ "sd"][ *col_itr] = mean_sd_collectorb.GetStandardDeviation();
      }

      // end
      return stats;
    }

    //! @brief initialize command line with all static flags and parameters
    //! @return Command object containing all the static flags and parameters from the FitInDensity class
    //! used to define its behavior
    util::ShPtr< command::Command> ProteinScore::InitializeCommand() const
    {
      // initialize a new command
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add all the flags and parameters
      sp_cmd->AddFlag( m_PDBFileFlag);

      sp_cmd->AddFlag( m_PDBFileListFlag);
      sp_cmd->AddFlag( m_DetailedFlag);
      sp_cmd->AddFlag( m_ScoreTableReadFlag);
      sp_cmd->AddFlag( m_ScoreTableWriteFlag);
      sp_cmd->AddFlag( fold::DefaultFlags::GetFlagNativeModel());
      sp_cmd->AddFlag( quality::Measures::GetFlagQualityMeasures());
      sp_cmd->AddFlag( biol::AtomTypes::GetFlagAtomTypes());
      sp_cmd->AddFlag( m_EnrichmentFlag);
      sp_cmd->AddFlag( m_WeightSetFlag);
      sp_cmd->AddFlag( m_IdealizeFlag);
      sp_cmd->AddFlag( score::Symmetry< assemble::ProteinModel>::GetFlagScoreSymmetry());
      sp_cmd->AddFlag( m_ContactFlag);

      sp_cmd->PushBack( fold::ProtocolRestraint::GetInstance().GetAllFlags());

      sp_cmd->AddFlag( fold::ProtocolEM::GetFlagScoreDensityAgreement());
      sp_cmd->AddFlag( fold::ProtocolEM::GetFlagBodyRestraint());
      sp_cmd->AddFlag( fold::ProtocolEM::GetFlagScoreDensityConnectivity());
      sp_cmd->AddFlag( m_DensityProfileAgreementFlag);

      sp_cmd->AddFlag( sspred::Methods::GetFlagReadSSPredictions());
      sp_cmd->AddFlag( fold::DefaultFlags::GetFlagReadSequenceDataPath());

      sp_cmd->AddFlag( score::AANeighborhoodExposurePrediction::GetFlagScoreExposure());

      //flag and three parameters adjusting the membrane
      sp_cmd->AddFlag( m_PDBTM_XML_flag);
      sp_cmd->AddFlag( biol::Membrane::GetFlagMembrane());

      // add the global flags
      sp_cmd->AddFlag( pdb::Handler::GetFlagHelixClasses());
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());
      sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());
      sp_cmd->AddFlag( pdb::Factory::GetFlagConvertToNaturalAAType());
      sp_cmd->AddFlag( pdb::Factory::GetFlagBiomolecule());
      sp_cmd->AddFlag( pdb::Factory::GetFlagDSSP());
      sp_cmd->AddFlag( m_DSSPNativeOnlyFlag);
      sp_cmd->AddFlag( quality::SuperimposeMeasures::GetFlagSuperimposeMeasure());

      sp_cmd->AddFlag( score::ProteinModelMembraneTopology::GetFlagExpectedTransmembraneHelicesPoolFile());

      sp_cmd->AddFlag( assemble::SSEPool::GetFlagPoolRead());
      sp_cmd->AddFlag( fold::DefaultFlags::GetFlagUseNativeSSEsAsPool());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief default constructor
    ProteinScore::ProteinScore() :
      m_PDBFileFlag
      (
        new command::FlagStatic
        (
          pdb::GetDefaultFileExtension(), "single pdb file",
          command::Parameter
          (
            "pdb_filename",
            "\tfilename for input pdb to be scored",
            command::ParameterCheckOr
            (
              command::ParameterCheckExtension( ".pdb"),
              command::ParameterCheckExtension( ".pdb.bz2"),
              command::ParameterCheckExtension( ".pdb.gz"),
              command::ParameterCheckExtension( ".ent.gz")
            ),
            ""
          )
        )
      ),
      m_PDBFileListFlag
      (
        new command::FlagStatic
        (
          "pdblist",
          "file with pdblist",
          command::Parameter
          (
            "pdb_listfile",
            "\tfilename of a list of pdbfiles to be scored",
            command::ParameterCheckExtension( ".ls"),
            ""
          )
        )
      ),
      m_ListRangeParam
      (
        new command::Parameter
        (
          "list_range",
          "range of elements to use from the list; can be used if multiple score proteins need to be started for that "
          "same list, where the entire list would take too long to score",
          math::Range< size_t>( 0, util::GetUndefined< size_t>()).GetString()
        )
      ),
      m_PDBTM_XML_flag
      (
        new command::FlagStatic
        (
          "pdbtm_xml",
          "indicate that listfile contains the pdbtm xml file location after each pdb file location. This will be used to define a transformation that is applied to the pdb and contains the membrane core thickness, that overwrites the core thickness given in the commandline. It also forces the generation of membrane score objects."
        )
      ),
      m_DetailedFlag( new command::FlagStatic( "detailed", "write details of scores")),
      m_ScoreTableReadFlag
      (
        new command::FlagStatic
        (
          "score_table_read",
          "read a score table instead of calculating scores",
          command::Parameter( "score_table_name", "filename where score table gets read from", "scores.table")
        )
      ),
      m_ScoreTableWriteFlag
      (
        new command::FlagStatic
        (
          "score_table_write",
          "write all score to a given table",
          command::Parameter( "score_table_name", "filename where score table gets written to", "scores.table")
        )
      ),
      m_EnrichmentFlag( new command::FlagStatic( "enrichment", "calculate enrichment for scores")),
      m_EnrichmentGoodModelsFractionParam
      (
        new command::Parameter
        (
          "fraction_good_models",
          "define fraction of good models vs total number models",
          command::ParameterCheckRanged< double>( 0.0, 1.0),
          "0.1"
        )
      ),
      m_EnrichmentGoodModelsCriteriaCutoffParam
      (
        new command::Parameter( "criteria_cutoff", "define criteria of good models", "5.0")
      ),
      m_EnrichmentNumberCrossValidationParam
      (
        new command::Parameter
        (
          "number_cross_validation",
          "to achieve the fraction either good or bad models get removed systematically. To vary that and later calculate a standard deviation, the enrichment calculation will be repeated n times",
          command::ParameterCheckRanged< size_t>( 1, 100),
          "10"
        )
      ),
      m_EnrichmentCriteriaNameParam( new command::Parameter( "criteria", "column name for criteria", "RMSD100")),
      m_EnrichmentSortOrderParam
      (
        new command::Parameter
        (
          "order",
          "a selection of order functions",
          command::ParameterCheckEnumerate< math::Comparisons< double> >(),
          math::Comparisons< double>::GetEnums().e_Less.GetName()
        )
      ),
      m_EnrichmentOutputTableFilenameParam
      (
        new command::Parameter
        (
          "enrichment_table_filename", "filename of table to be outputted containign the enrichments", "enrichments.tbl"
        )
      ),
      m_EnrichmentPlotFilePrefixParam
      (
        new command::Parameter( "plot_file_prefix", "if given, write roc plots to files with this prefix", "")
      ),
      m_WeightSetFlag
      (
        new command::FlagStatic
        (
          "weight_set",
          "file containing storage::Table< double> with weight set",
          command::Parameter( "weight_set_file_name", "filename of file containing weight set", "")
        )
      ),
      m_IdealizeFlag( new command::FlagStatic( "idealize", "idealize pdb before scoring")),
      m_ContactFlag( new command::FlagStatic( "contact", "\tcalculate contact statistics and contact recovery values")),
      m_DSSPNativeOnlyFlag
      (
        new command::FlagStatic
        (
          "dssp_native",
          "Use DSSP only to assign sses from the native structure. This is appropriate if using models generated by "
          "bcl which may have some strands too far away from on another or too short (w/o loops) for DSSP to classify "
          "them as strands"
        )
      ),
      m_DensityProfileAgreementFlag
      (
        new command::FlagStaticAndDynamic
        (
          "density_profile_agreement",
          "score density profiles agreement - native is used as body restraint",
          command::Parameter
          (
            "profile",
            "the profiles to be used for scoring",
            command::ParameterCheckSerializable( score::DensityProfileSSEAgreement::Profile1DTypeEnum())
          ),
          0,
          score::DensityProfileSSEAgreement::s_NumberProfile1DTypes
        )
      ),
      m_MRCMapParam
      (
        new command::Parameter
        (
          "mrc_filename",
          "\tfilename for evaluating fit of protein",
          command::ParameterCheckFileExistence(),
          ""
        )
      ),
      m_MRCResolutionParam
      (
        new command::Parameter
        (
          "density map resolution",
          "\tresolution of given electron density [A] - will be used to simulate density and calculate correlation",
          command::ParameterCheckRanged< double>( 0.0, 100.0),
          "0.0"
        )
      ),
      m_SimulatorParameter
      (
        new command::Parameter
        (
          "simulator",
          "the simulator to be used to generate density maps from a atom structure",
          command::ParameterCheckEnumerate< density::Simulators>(),
          density::GetSimulators().e_Gaussian.GetName()
        )
      ),
      m_Membrane( new biol::Membrane()),
      m_ScoringTable(),
      m_Weights(),
      m_TemplateModel(),
      m_ScoreFunction(),
      m_AtomTypes(),
      m_Qualities(),
      m_ProteinModelData( new assemble::ProteinModelData())
    {
      // add parameters to flags
      m_PDBFileListFlag->PushBack( m_ListRangeParam);
      m_EnrichmentFlag->PushBack( m_EnrichmentGoodModelsFractionParam);
      m_EnrichmentFlag->PushBack( m_EnrichmentGoodModelsCriteriaCutoffParam);
      m_EnrichmentFlag->PushBack( m_EnrichmentNumberCrossValidationParam);
      m_EnrichmentFlag->PushBack( m_EnrichmentCriteriaNameParam);
      m_EnrichmentFlag->PushBack( m_EnrichmentSortOrderParam);
      m_EnrichmentFlag->PushBack( m_EnrichmentOutputTableFilenameParam);
      m_EnrichmentFlag->PushBack( m_EnrichmentPlotFilePrefixParam);
      m_DensityProfileAgreementFlag->PushBack( m_MRCMapParam);
      m_DensityProfileAgreementFlag->PushBack( m_MRCResolutionParam);
      m_DensityProfileAgreementFlag->PushBack( m_SimulatorParameter);
    }

  } // namespace app
} // namespace bcl
